# Spatial Basis Density Field with Boundary Attenuation

This note describes the experimental density model implemented by
`spatial_basis_density_field()`.

The goal is to estimate a continuous transcript density field after
transcript-level cell-type posterior probabilities have already been inferred,
for example by `spatial_basis_segmentation()`.

The default mode fits one shared total transcript density field. Cell-type
specific densities are then computed by multiplying this total density by the
input posterior probabilities. This simpler model is often preferable when the
main goal is to find transcript-poor boundaries shared across cell types.

## Observed Data

For each transcript `i`, we use:

- location `x_i`
- posterior cell-type weights `q_ik`

where:

$$
q_{ik} = P(\text{cell type} = k \mid \text{transcript } i)
$$

In the default shared-density mode, the posterior weights are not used to fit
separate density fields. Instead, all transcripts contribute to one total
transcript point process, and the posteriors are used afterward to decompose
the total density by cell type.

The older cell-type-specific mode treats the posterior weights as fractional
observations of cell-type-specific transcript point processes.

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
\lambda_{ia} \ge 0
$$

In 2D, each point has three active barycentric coordinates. In 3D, each point
has four.

## Density Model

The default shared-density model fits one total transcript density:

$$
\rho(x)
= \exp(\eta(x))
  \left[
    f + (1 - f)\,\sigma(-g(x))
  \right]
$$

Here:

- `eta(x)` is a broad baseline log-density field
- `g(x)` is a boundary logit field
- `f`, called `floor` or `density_floor` in the implementation, is a fixed
  minimum attenuation multiplier

The attenuation term:

$$
A(x) = f + (1 - f)\,\sigma(-g(x))
$$

is near `1` when `g(x)` is very negative, and near `density_floor` when `g(x)` is
very positive.

Thus high positive values of `g(x)` correspond to boundary-like regions where
the total transcript density is attenuated.

Cell-type-specific density is then:

$$
\rho_k(x_i) = \rho(x_i) q_{ik}
$$

where $q_{ik}$ is the supplied transcript-level posterior probability. This is
the value returned in the `density` matrix when `density_mode = "shared"`.
The shared total density is returned separately as `total_density`.

The optional cell-type-specific mode can still be requested with:

```r
density_mode = "cell_type"
```

In that mode, each cell type gets its own density field:

$$
\rho_k(x)
= \exp(\eta_k(x))
  \left[
    f_k + (1 - f_k)\,\sigma(-g_k(x))
  \right].
$$

The floor prevents boundaries from driving density exactly to zero. For example:

```r
density_floor = 0.05
```

means a strong boundary can reduce density to 5 percent of the baseline
density, but not below that.

## Baseline Density Field

The baseline log-density is vertex-linear:

$$
\eta(x_i) = \sum_a \lambda_{ia} w_a
$$

where `w_a` is the fitted baseline log-density coefficient for basis vertex
`a`. In `density_mode = "cell_type"`, this becomes $w_{a,k}$.

This is the smooth background density field before boundary attenuation.

## Boundary Logit Field

The boundary logit field combines vertex-linear terms with edge-quadratic
terms:

$$
g(x_i)
= \sum_a \lambda_{ia} s_a
  + \sum_{a<b} \lambda_{ia}\lambda_{ib}s_{ab}
$$

The vertex terms:

$$
\sum_a \lambda_{ia} s_a
$$

define a continuous piecewise-linear boundary field.

The edge terms:

$$
\lambda_{ia}\lambda_{ib}s_{ab}
$$

allow the boundary logit to bend or form a localized interior high-probability
region within a simplex. These are the terms that let the model represent a
thin low-density band between high-density regions.

The fitted edge coefficient `s_ab` is shared globally for the active basis
pair `(a, b)`. Sharing these coefficients across adjacent simplexes helps
preserve continuity of the boundary contribution.

In `density_mode = "cell_type"`, the corresponding coefficients are
cell-type-specific: $s_{a,k}$ and $s_{ab,k}$.

## Interpretation of `g`

The model uses:

$$
A(x) = f + (1 - f)\,\sigma(-g(x))
$$

so:

$$
\begin{aligned}
g(x) &\ll 0 &&\Rightarrow A(x) \approx 1, \\
g(x) &\gg 0 &&\Rightarrow A(x) \approx f, \\
g(x) &= 0   &&\Rightarrow A(x) = \frac{1 + f}{2}.
\end{aligned}
$$

In other words, `sigmoid(g(x))` can be interpreted as a boundary probability
or boundary propensity, while `sigmoid(-g(x))` is the density-preservation
factor before applying the floor.

## Likelihood

In shared-density mode, all transcripts are modeled as one inhomogeneous
Poisson point process with intensity:

$$
\rho(x)
$$

The weighted negative log-likelihood is:

$$
\mathcal{L}_{\text{data}}
= \int_\Omega \rho(x)\,dx
  - \sum_i \log \rho(x_i)
$$

up to constants that do not depend on the fitted density field.

In `density_mode = "cell_type"`, the older posterior-weighted likelihood is
used:

$$
\mathcal{L}_{\text{data}}
= \sum_k \int_\Omega \rho_k(x)\,dx
  - \sum_i \sum_k q_{ik} \log \rho_k(x_i).
$$

The implementation divides the objective and gradient by the total posterior
mass:

$$
Q = n
$$

in shared mode, or $Q = \sum_i \sum_k q_{ik}$ in cell-type mode, so the scale
is roughly per transcript.

## Quadrature Approximation

The integral:

$$
\int_\Omega \rho_k(x)\,dx
$$

is approximated numerically.

In shared mode, the corresponding integral is $\int_\Omega \rho(x)\,dx$.

The default implementation uses occupied-simplex quadrature. It first finds the
spatial basis simplex containing each transcript, keeps the unique occupied
simplexes, and places quadrature points inside each one.

For `basis = "2d"`, `quadrature_subdivision = m` splits each occupied triangle
into `m^2` equal-area subtriangles and places one quadrature point at each
subtriangle centroid. Each quadrature point has weight:

$$
v_j = \frac{\operatorname{area}(\Delta)}{m^2}.
$$

For `basis = "3d"`, `quadrature_subdivision = m` currently uses `m^3`
simplex-local Duffy-midpoint quadrature points per occupied tetrahedron. The
weights are scaled so they sum to the tetrahedron volume.

If quadrature points $z_j$ have weights $v_j$, then:

$$
\int_\Omega \rho_k(x)\,dx
\approx
\sum_j v_j \rho_k(z_j)
$$

This simplex-local quadrature is usually a better match to the model than a
rectangular grid because the density field is defined by barycentric
coordinates inside simplexes. It also reduces aliasing artifacts where the
optimizer can create very sharp density peaks or dips between sparse
rectangular grid points.

The older rectangular-grid quadrature can still be requested with:

```r
quadrature_method = "grid"
```

The user controls the approximate number of quadrature points with:

```r
quadrature_n
```

and the fractional bounding-box expansion with:

```r
quadrature_expansion
```

These grid-specific arguments are only used when `quadrature_method = "grid"`.

### Adaptive Quadrature

An adaptive simplex grid would be a natural future extension. The basic idea
would be:

1. Fit the model with a coarse simplex quadrature.
2. Evaluate the fitted field inside each occupied simplex.
3. Refine simplexes where the fitted field has large variation, large gradient
   magnitude, high curvature, or extreme transcript-versus-quadrature density
   mismatch.
4. Refit using the refined quadrature rule.

This would target quadrature effort near thin boundaries or sharp peaks,
instead of spending the same number of points in every occupied simplex. The
important implementation detail is that the quadrature rule should be fixed
during each optimization run; changing quadrature points continuously during
L-BFGS would make the objective unstable.

## Regularization

The model has several regularization terms.

### Baseline Density Smoothness

The baseline log-density coefficients are smoothed over neighboring lattice
basis points:

$$
R_\eta
= \lambda_\eta
  \sum_{(a,b) \in E}
  \frac{1}{2}
  \left(
    \frac{w_a - w_b}{d_{ab}}
  \right)^2
$$

This discourages sharp baseline density variation.

### Boundary Vertex Smoothness

The vertex boundary logits are also smoothed:

$$
R_{s,\text{smooth}}
= \lambda_{s,\text{smooth}}
  \sum_{(a,b) \in E}
  \frac{1}{2}
  \left(
    \frac{s_a - s_b}{d_{ab}}
  \right)^2
$$

This discourages highly oscillatory boundary fields.

### Boundary Prior

The vertex boundary logits are shrunk toward:

```r
s_prior_mean
```

with strength:

```r
lambda_s_prior
```

A negative value such as:

```r
s_prior_mean = -4
```

favors no boundary attenuation unless the data support positive `g(x)`.

The corresponding penalty is:

$$
R_{s,\text{prior}}
= \lambda_{s,\text{prior}}
  \sum_a
  \frac{1}{2}
  \left(s_a - \mu_s\right)^2
$$

where $\mu_s$ is `s_prior_mean`.

### Edge-Quadratic Shrinkage

The edge-quadratic boundary terms are shrunk toward zero:

$$
R_{\text{edge}}
= \lambda_{\text{edge}}
  \sum_{(a,b)}
  \frac{1}{2}s_{ab}^2
$$

This is important because these terms can create localized dips and should not
be allowed to create spurious holes everywhere.

The fitted objective is the posterior-mass-scaled data term plus these
regularizers:

$$
\mathcal{J}
= \frac{1}{Q}\mathcal{L}_{\text{data}}
  + R_\eta
  + R_{s,\text{smooth}}
  + R_{s,\text{prior}}
  + R_{\text{edge}}.
$$

## Parameters

The fitted parameter vector contains:

```text
w_a        baseline log-density vertex coefficients
s_a        boundary-logit vertex coefficients
s_ab       boundary-logit edge-quadratic coefficients
```

In `density_mode = "cell_type"`, these become cell-type-specific coefficients:
`w_a,k`, `s_a,k`, and `s_ab,k`.

The returned object exposes these as:

```r
eta_basis
boundary_node_basis
boundary_edge_basis
boundary_edge_pairs
```

Transcript-level predictions are returned as:

```r
density
total_density
eta
boundary_logit
attenuation
posterior
```

## Example

```r
density_fit <- spatial_basis_density_field(
    transcripts_df = fit$transcripts_df,
    posterior = fit$marginals,
    basis = "2d",
    s = 2,
    density_mode = "shared",
    density_floor = 0.05,
    quadrature_method = "simplex",
    quadrature_subdivision = 4,
    lambda_eta = 0.1,
    lambda_s_smooth = 0.1,
    lambda_s_prior = 0.1,
    s_prior_mean = -4,
    lambda_edge = 0.25,
    maxit = 100
)
```

## Current Limitations

This is a prototype. Important limitations include:

- The integral is approximated over a rectangular bounding box, not the true
  tissue support when `quadrature_method = "grid"`. The default simplex
  quadrature only integrates over occupied simplexes.
- The floor is fixed, not learned.
- Boundary attenuation can be partially confounded with the baseline density
  field.
- Edge-quadratic terms increase flexibility but also increase the need for
  regularization.
- In shared-density mode, cell-type density is only as spatially resolved as
  the supplied posterior probabilities.
- The model estimates density fields, not cell partitions directly.

The intended next step is to inspect fitted density and boundary fields and
decide whether they provide useful information for downstream transcript
partitioning.
