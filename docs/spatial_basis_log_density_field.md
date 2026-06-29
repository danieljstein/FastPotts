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

Each transcript or quadrature point has a small set of active basis vertices:

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

## Simplex-Local Quadrature

The integral:

$$
\int_\Omega \exp(\eta(x))\,dx
$$

is approximated using occupied-simplex quadrature. The implementation first
finds the spatial basis simplex containing each transcript, keeps the unique
occupied simplexes, and places quadrature points inside each one.

For `basis = "2d"`, `quadrature_subdivision = m` splits each occupied triangle
into `m^2` equal-area subtriangles and places one quadrature point at each
subtriangle centroid. Each quadrature point has weight:

$$
v_j = \frac{\operatorname{area}(\Delta)}{m^2}.
$$

For `basis = "3d"`, `quadrature_subdivision = m` currently uses `m^3`
simplex-local Duffy-midpoint quadrature points per occupied tetrahedron. The
weights are scaled so they sum to the tetrahedron volume.

If quadrature points `z_j` have weights `v_j`, then:

$$
\int_\Omega \exp(\eta(x))\,dx
\approx
\sum_j v_j \exp(\eta(z_j)).
$$

This keeps the integral approximation aligned with the barycentric basis used
by the density model.

## Regularization

The model regularizes neighboring basis-point log-density coefficients. For a
neighboring basis edge `(a, b)`, define the spatial slope:

$$
r_{ab} = \frac{w_a - w_b}{d_{ab}},
$$

where `d_ab` is the physical distance between neighboring basis points.

The regularization term is:

$$
R(w)
= \lambda
  \frac{1}{|E|}
  \sum_{(a,b) \in E}
  \psi(r_{ab}).
$$

The fitted objective is:

$$
\mathcal{J}
= \frac{1}{n}\mathcal{L}_{\text{data}} + R(w).
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

## Example

```r
log_density_fit <- spatial_basis_log_density_field(
    transcripts_df = fit$transcripts_df,
    posterior = fit$marginals,
    basis = "2d",
    s = 2,
    quadrature_subdivision = 4,
    lambda = 0.1,
    regularization = "huber",
    delta = 1,
    sigma = 1,
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

- The integral is only over occupied simplexes, not a full tissue mask.
- Very fine meshes may still require stronger regularization to avoid
  overfitting.
- Cell-type-specific density is only as spatially resolved as the supplied
  posterior probabilities.
- The model estimates density fields, not cell partitions directly.
