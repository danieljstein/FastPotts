# Spatial Basis Density Field with Boundary Attenuation

This note describes the experimental density model implemented by
`spatial_basis_density_field()`.

The goal is to estimate a continuous transcript density field for each cell
type after transcript-level cell-type posterior probabilities have already been
inferred, for example by `spatial_basis_segmentation()`.

## Observed Data

For each transcript `i`, we use:

- location `x_i`
- posterior cell-type weights `q_ik`

where:

$$
q_{ik} = P(\text{cell type} = k \mid \text{transcript } i)
$$

The model treats the posterior weights as fractional observations of cell-type
specific transcript point processes.

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

For cell type `k`, the modeled transcript density is:

$$
\rho_k(x)
= \exp(\eta_k(x))
  \left[
    f_k + (1 - f_k)\,\sigma(-g_k(x))
  \right]
$$

Here:

- `eta_k(x)` is a broad baseline log-density field
- `g_k(x)` is a boundary logit field
- `f_k`, called `floor_k` in the implementation, is a fixed
  cell-type-specific minimum attenuation multiplier

The attenuation term:

$$
A_k(x) = f_k + (1 - f_k)\,\sigma(-g_k(x))
$$

is near `1` when `g_k(x)` is very negative, and near `floor_k` when `g_k(x)` is
very positive.

Thus high positive values of `g_k(x)` correspond to boundary-like regions where
the density of type `k` is attenuated.

The floor prevents boundaries from driving density exactly to zero. For example:

```r
density_floor = 0.05
```

means a strong boundary can reduce density to 5 percent of the baseline
density, but not below that.

## Baseline Density Field

The baseline log-density is vertex-linear:

$$
\eta_k(x_i) = \sum_a \lambda_{ia} w_{a,k}
$$

where `w_{a,k}` is the fitted baseline log-density coefficient for basis vertex
`a` and cell type `k`.

This is the smooth background density field before boundary attenuation.

## Boundary Logit Field

The boundary logit field combines vertex-linear terms with edge-quadratic
terms:

$$
g_k(x_i)
= \sum_a \lambda_{ia} s_{a,k}
  + \sum_{a<b} \lambda_{ia}\lambda_{ib}s_{ab,k}
$$

The vertex terms:

$$
\sum_a \lambda_{ia} s_{a,k}
$$

define a continuous piecewise-linear boundary field.

The edge terms:

$$
\lambda_{ia}\lambda_{ib}s_{ab,k}
$$

allow the boundary logit to bend or form a localized interior high-probability
region within a simplex. These are the terms that let the model represent a
thin low-density band between high-density regions.

The fitted edge coefficient `s_ab,k` is shared globally for the active basis
pair `(a, b)`. Sharing these coefficients across adjacent simplexes helps
preserve continuity of the boundary contribution.

## Interpretation of `g`

The model uses:

$$
A_k(x) = f_k + (1 - f_k)\,\sigma(-g_k(x))
$$

so:

$$
\begin{aligned}
g_k(x) &\ll 0 &&\Rightarrow A_k(x) \approx 1, \\
g_k(x) &\gg 0 &&\Rightarrow A_k(x) \approx f_k, \\
g_k(x) &= 0   &&\Rightarrow A_k(x) = \frac{1 + f_k}{2}.
\end{aligned}
$$

In other words, `sigmoid(g_k(x))` can be interpreted as a boundary probability
or boundary propensity, while `sigmoid(-g_k(x))` is the density-preservation
factor before applying the floor.

## Likelihood

For each cell type `k`, transcripts are modeled as an inhomogeneous Poisson
point process with intensity:

$$
\rho_k(x)
$$

Because transcript cell types are uncertain, each transcript contributes with
posterior weight `q_ik`.

The weighted negative log-likelihood is:

$$
\mathcal{L}_{\text{data}}
= \sum_k \int_\Omega \rho_k(x)\,dx
  - \sum_i \sum_k q_{ik} \log \rho_k(x_i)
$$

up to constants that do not depend on the fitted density field.

The implementation divides the objective and gradient by the total posterior
mass:

$$
Q = \sum_i \sum_k q_{ik}
$$

so the scale is roughly per transcript.

## Quadrature Approximation

The integral:

$$
\int_\Omega \rho_k(x)\,dx
$$

is approximated numerically.

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
  \sum_k
  \frac{1}{2}
  \left(
    \frac{w_{a,k} - w_{b,k}}{d_{ab}}
  \right)^2
$$

This discourages sharp baseline density variation.

### Boundary Vertex Smoothness

The vertex boundary logits are also smoothed:

$$
R_{s,\text{smooth}}
= \lambda_{s,\text{smooth}}
  \sum_{(a,b) \in E}
  \sum_k
  \frac{1}{2}
  \left(
    \frac{s_{a,k} - s_{b,k}}{d_{ab}}
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

favors no boundary attenuation unless the data support positive `g_k(x)`.

The corresponding penalty is:

$$
R_{s,\text{prior}}
= \lambda_{s,\text{prior}}
  \sum_a \sum_k
  \frac{1}{2}
  \left(s_{a,k} - \mu_s\right)^2
$$

where $\mu_s$ is `s_prior_mean`.

### Edge-Quadratic Shrinkage

The edge-quadratic boundary terms are shrunk toward zero:

$$
R_{\text{edge}}
= \lambda_{\text{edge}}
  \sum_{(a,b)} \sum_k
  \frac{1}{2}s_{ab,k}^2
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

The fitted parameter vector contains, for every cell type:

```text
w_a,k      baseline log-density vertex coefficients
s_a,k      boundary-logit vertex coefficients
s_ab,k     boundary-logit edge-quadratic coefficients
```

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
eta
boundary_logit
attenuation
```

## Example

```r
density_fit <- spatial_basis_density_field(
    transcripts_df = fit$transcripts_df,
    posterior = fit$marginals,
    basis = "2d",
    s = 2,
    density_floor = 0.05,
    quadrature_n = 10000,
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
  tissue support.
- The floor is fixed, not learned.
- Boundary attenuation can be partially confounded with the baseline density
  field.
- Edge-quadratic terms increase flexibility but also increase the need for
  regularization.
- The model estimates density fields, not cell partitions directly.

The intended next step is to inspect fitted density and boundary fields and
decide whether they provide useful information for downstream transcript
partitioning.
