# Spatial Basis Log-Density Field

This note describes the simplified density model implemented by
`spatial_basis_log_density_field()`.

The goal is to estimate a continuous total transcript density field without
separate boundary probability or attenuation terms. This model is intended as a
cleaner baseline for finding transcript-poor regions using a fine spatial mesh
and edge-preserving regularization.

If transcript-level cell-type posterior probabilities are supplied, they are
not used to fit separate density fields. Instead, the model fits one total
transcript density field and then decomposes that density by cell type:

$$
\rho_k(x_i) = \rho(x_i) q_{ik}.
$$

For large datasets, this derived transcript-by-cell-type matrix can be much
larger than the fitted total density vector. If only the total density is
needed, use:

```r
return_density = FALSE
```

This still allows `posterior` to be supplied and normalized, but skips returning
the dense `density` matrix.

## Observed Data

For each transcript `i`, the model uses its spatial location:

$$
x_i.
$$

Optionally, the user can also provide posterior cell-type weights:

$$
q_{ik} = P(\text{cell type} = k \mid \text{transcript } i).
$$

These posterior weights are only used after fitting, to compute
cell-type-specific transcript densities from the fitted total density.

## Spatial Basis

The model uses the same barycentric lattice basis as
`spatial_basis_segmentation()`:

- `basis = "2d"`: triangular lattice
- `basis = "3d"`: body-centered cubic lattice

Each transcript, prediction point, or integration point has a small set of
active basis vertices:

$$
\lambda_{ia}
$$

with:

$$
\sum_a \lambda_{ia} = 1,
\qquad
\lambda_{ia} \ge 0.
$$

In 2D, each point has three active barycentric coordinates. In 3D, each point
has four.

## Density Model

The fitted total transcript density is:

$$
\rho(x) = \exp(\eta(x)).
$$

The log-density field is vertex-linear within each simplex:

$$
\eta(x_i) = \sum_a \lambda_{ia} w_a,
$$

where `w_a` is the fitted log-density coefficient for basis vertex `a`.

There are no boundary logits, attenuation floors, or edge-quadratic terms in
this model. Low-density boundaries must be represented directly as lower values
of `eta`.

## Likelihood

Transcripts are modeled as an inhomogeneous Poisson point process with
intensity:

$$
\rho(x).
$$

The negative log-likelihood is:

$$
\mathcal{L}_{\text{data}}
= \int_\Omega \rho(x)\,dx
  - \sum_i \log \rho(x_i)
$$

up to constants that do not depend on the fitted density field.

Substituting the log-linear model gives:

$$
\mathcal{L}_{\text{data}}
= \int_\Omega \exp(\eta(x))\,dx
  - \sum_i \eta(x_i).
$$

The implementation divides the objective and gradient by the number of
transcripts, so the data term is roughly on a per-transcript scale.

## Integration Domain and Simplex Integration

The integral:

$$
\int_\Omega \exp(\eta(x))\,dx
$$

is computed over a simplex domain built from the spatial basis. By default,
the implementation finds the basis simplex containing each transcript and keeps
the unique occupied simplexes.

The domain can also be expanded before integration:

```r
domain_expansion_steps = 2L
domain_expansion_axes = "xy"
```

`domain_expansion_steps = 0L` keeps the occupied-simplex-only behavior. Larger
values add basis vertices within that many basis-graph steps of the occupied
vertices, then include all simplexes whose vertices are in the expanded basis
set. In 3D, `domain_expansion_axes = "xy"` expands only along graph edges with
zero z displacement, which is useful for thin samples where empty space should
be filled laterally but not through z.

For each simplex $\Delta$, the contribution to the normalizing integral is:

$$
I_\Delta(w)
=
\int_\Delta \exp\left(\sum_a \lambda_a(x) w_a\right)\,dx.
$$

The default integration method is:

```r
integration_method = "subdivision_linear"
```

This uses a positive piecewise-linear interpolation approximation over the
barycentric subdivision of each simplex:

$$
I_\Delta(w)
\approx
|\Delta|
\sum_{q \in Q}
\omega_q
\exp\left(\sum_a \lambda_{qa} w_a\right).
$$

The corresponding gradient contribution is:

$$
\frac{\partial I_\Delta}{\partial w_a}
\approx
|\Delta|
\sum_{q \in Q}
\omega_q
\lambda_{qa}
\exp\left(\sum_b \lambda_{qb} w_b\right).
$$

For 2D triangles, the subdivision-linear rule uses these keypoint weights:

- original vertices: weight $1/9$ each
- edge midpoints: weight $1/9$ each
- simplex centroid: weight $1/3$

For 3D tetrahedra, it uses:

- original vertices: weight $1/16$ each
- edge midpoints: weight $1/24$ each
- face centroids: weight $1/16$ each
- simplex centroid: weight $1/4$

These weights are normalized within each simplex, so they sum to one before
multiplication by the simplex volume. The approximation is fast, positive, and
compatible with the gradient used by L-BFGS. Because $\exp(\eta)$ is convex in
the log-linear field, this piecewise-linear approximation tends to overestimate
the exact simplex integral.

For higher accuracy, use:

```r
integration_method = "analytic"
```

The analytic method computes the exact log-linear simplex integral with stable
Taylor and Gauss fallback paths for near-identical vertex values. It is usually
slower than `subdivision_linear`, especially for expanded domains with many
empty simplexes.

The compatibility arguments:

```r
quadrature_subdivision
store_quadrature_coords
```

are retained for older code paths but are not used by the current
`spatial_basis_log_density_field()` objective.

## Regularization

The model regularizes neighboring basis-point log-density coefficients. For a
neighboring basis edge `(a, b)`, define the spatial slope:

$$
r_{ab} = \frac{w_a - w_b}{d_{ab}},
$$

where `d_ab` is the physical distance between neighboring basis points.

The regularization term is:

$$
R_{\nabla}(w)
= \lambda
  \frac{1}{|E|}
  \sum_{(a,b) \in E}
  \psi(r_{ab}).
$$

The fitted objective is:

$$
\mathcal{J}
= \frac{1}{n}\mathcal{L}_{\text{data}}
  + R_{\nabla}(w)
  + R_{\Delta}(w).
$$

### Quadratic

With `regularization = "quadratic"`:

$$
\psi(r) = \frac{1}{2}r^2.
$$

This is the smoothest and most stable option, but it tends to blur sharp
density changes.

### Huber

With `regularization = "huber"`:

$$
\psi(r) =
\begin{cases}
\frac{1}{2}r^2, & |r| \le \delta, \\
\delta(|r| - \frac{1}{2}\delta), & |r| > \delta.
\end{cases}
$$

This behaves quadratically for small slopes and linearly for large slopes. It
is a useful first edge-preserving option.

### Bounded

With `regularization = "bounded"`:

$$
\psi(r)
= \sigma^2
  \left[
    1 - \frac{1}{1 + (r/\sigma)^2}
  \right].
$$

This has quadratic small-slope behavior but saturates for large slopes. It can
preserve stronger edges, but the objective is more nonconvex.

### Graph Laplacian Curvature Penalty

The optional `lambda_laplacian` penalty acts on the graph Laplacian of the
basis-point log-density field. For each basis point `a`, define a
distance-weighted neighbor average:

$$
\bar{w}_a =
\frac{\sum_{b \in N(a)} \omega_{ab} w_b}
     {\sum_{b \in N(a)} \omega_{ab}},
\qquad
\omega_{ab} = \frac{1}{d_{ab}^2}.
$$

The discrete Laplacian residual is:

$$
(\Delta w)_a = w_a - \bar{w}_a.
$$

The penalty is:

$$
R_{\Delta}(w)
= \lambda_{\Delta}
  \frac{1}{2M}
  \sum_a
  \left[(\Delta w)_a\right]^2,
$$

where `M` is the number of basis points.

This penalty suppresses isolated speckles and jagged curvature by encouraging
each basis coefficient to agree with a local linear/harmonic continuation from
its neighbors. It is different from the first-difference penalty: the
first-difference penalty suppresses slopes, while the Laplacian penalty
suppresses changes in slope.

The two penalties can be combined. A useful starting point is:

```r
regularization = "bounded"
lambda = 0.1
lambda_laplacian = 0.1
```

## Parameters

The fitted parameter vector contains one coefficient per basis point:

```text
w_a      log-density vertex coefficient
```

The returned object exposes these as:

```r
eta_basis
```

Transcript-level predictions are returned as:

```r
total_density
eta
```

If `posterior` is supplied, the returned object also contains:

```r
density
posterior
```

where:

$$
\rho_{ik} = \rho_i q_{ik}.
$$

The returned object also includes the basis lattice and integration domain:

```r
basis_lattice
basis_points
domain
```

For the current implementation, `domain$volume` is the scalar area or volume of
one simplex in the regular lattice, and `domain$n_simplex` is the number of
simplexes included in the integration domain. The effective options are stored
in:

```r
parameters
```

## Example

```r
log_density_fit <- spatial_basis_log_density_field(
    transcripts_df = fit$transcripts_df,
    posterior = fit$marginals,
    basis = "3d",
    s = 2,
    integration_method = "subdivision_linear",
    domain_expansion_steps = 2L,
    domain_expansion_axes = "xy",
    lambda = 1,
    regularization = "bounded",
    delta = 1,
    sigma = 1,
    lambda_laplacian = 100,
    maxit = 100
)
```

## When to Use This Model

This model is a good baseline when:

- the boundary-attenuation model is difficult to tune
- boundary logits are hard to interpret
- a fine mesh can directly represent transcript-poor bands
- the main goal is a stable total transcript density field

Compared with `spatial_basis_density_field()`, this model has fewer parameters
and less identifiability ambiguity. It cannot explicitly separate broad
baseline density from a multiplicative boundary attenuation term, but that
simplicity is often helpful.

## Current Limitations

Important limitations include:

- By default, the integral is over occupied simplexes. Domain expansion can add
  nearby empty simplexes, but it is still a graph-dilation heuristic rather than
  a full tissue mask.
- The default `subdivision_linear` integration is approximate. It is much
  faster than the analytic integral, but it tends to overestimate the
  normalizing integral for convex log-linear fields.
- Very fine meshes may still require stronger regularization to avoid
  overfitting.
- The graph Laplacian penalty suppresses speckles, but it does not explicitly
  distinguish blob-like curvature from line-like or sheet-like boundaries.
- Cell-type-specific density is only as spatially resolved as the supplied
  posterior probabilities.
- The model estimates density fields, not cell partitions directly.
