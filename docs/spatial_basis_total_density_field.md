# Spatial Basis Total Density Field

This note describes `spatial_basis_total_density_field()`.

The function fits a continuous total transcript density field on a spatial
basis that is aligned to an existing `spatial_basis_segmentation()` result. It
does not use posterior probabilities, cell-type labels, or gene signatures. Its
purpose is to estimate where transcripts are dense or sparse in space, while
reusing the segmentation basis as the tissue-supported domain.

This is intended as a cleaner input for basis-level cell partitioning. A later
watershed step can combine the fitted total density with the spatial prior from
`spatial_basis_segmentation()`:

$$
\rho_k(x) = \rho(x) p_k(x),
$$

where $\rho(x)$ is the total transcript density and $p_k(x)$ is the fitted
spatial prior for cell type $k$.

## Inputs

The main input is a fitted segmentation object:

```r
density_fit <- spatial_basis_total_density_field(
    segmentation_fit = fit,
    basis_subdivision = 2,
    quadrature_subdivision = 1,
    lambda = 1,
    lambda_laplacian = 1
)
```

The function uses:

- `fit$basis_lattice`
- `fit$basis_points`
- `fit$basis_edges`
- `fit$transcripts_df`
- `fit$parameters$basis`, `fit$parameters$s`, and `fit$parameters$origin` when available

For older segmentation objects without `fit$parameters`, the function falls
back to inferring `basis`, `s`, and `origin` from the basis geometry.

## Density Model

The fitted total transcript density is:

$$
\rho(x) = \exp(\eta(x)).
$$

The log-density field is piecewise linear on the density basis:

$$
\eta(x) = \sum_a \lambda_a(x) w_a,
$$

where:

- $\lambda_a(x)$ is the barycentric weight for density-basis vertex $a$
- $w_a$ is the fitted log-density coefficient at density-basis vertex $a$

The model ignores cell types. Each observed transcript contributes one point to
the total inhomogeneous Poisson point-process likelihood.

## Likelihood

The data term is the negative log-likelihood of an inhomogeneous Poisson point
process:

$$
\mathcal{L}_{\text{data}}
= \int_\Omega \rho(x)\,dx
  - \sum_i \log \rho(x_i),
$$

or equivalently:

$$
\mathcal{L}_{\text{data}}
= \int_\Omega \exp(\eta(x))\,dx
  - \sum_i \eta(x_i).
$$

The implementation divides the objective and gradient by the number of
transcripts, so the data term is on an approximate per-transcript scale.

## Parent Basis And Density Basis

There are two spatial bases:

1. The parent segmentation basis from `spatial_basis_segmentation()`
2. The density basis used to fit $\rho(x)$

The parent basis defines the integration domain $\Omega$. More concretely, the
function finds the parent simplexes that contain at least one transcript and
integrates over those parent simplexes.

The density basis may be the same as the parent basis, or a finer aligned
lattice controlled by:

```r
basis_subdivision
```

If the segmentation mesh size is $s$, the density mesh size is:

$$
s_{\text{density}} = \frac{s}{\text{basis\_subdivision}}.
$$

When `basis_subdivision = 1`, the density basis is the same basis as the
segmentation fit. In this case the returned `basis_points` and `basis_lattice`
line up row-for-row with `fit$basis_points` and `fit$basis_lattice`.

When `basis_subdivision > 1`, the density field is fit on a finer lattice. The
returned `basis_points` and `basis_lattice` then correspond to the density
basis, while `parent_basis_points` and `parent_basis_lattice` store the original
segmentation basis.

## Why Parent-Domain Quadrature Matters

The integral term:

$$
\int_\Omega \exp(\eta(x))\,dx
$$

is what makes absence of transcripts informative. If a region is inside the
tissue-supported parent basis but has no transcripts, it still contributes to
the integral. The optimizer can lower density there to reduce the integral
without sacrificing transcript likelihood.

This is why `spatial_basis_total_density_field()` integrates over occupied
parent simplexes, not only over density-basis simplexes that contain
transcripts.

## Basis Subdivision Vs Quadrature Subdivision

`basis_subdivision` changes the fitted model resolution.

`quadrature_subdivision` changes only the numerical integration resolution.

The effective quadrature subdivision inside each occupied parent simplex is:

$$
m_{\text{effective}}
= \text{basis\_subdivision}
  \times \text{quadrature\_subdivision}.
$$

For example:

```r
basis_subdivision = 3
quadrature_subdivision = 1
```

fits density on a $3\times$ finer basis and uses one quadrature layer per
density subdivision. For 2D, this gives approximately $m_{\text{effective}}^2$
quadrature points per occupied parent triangle. For 3D, this gives
approximately $m_{\text{effective}}^3$ quadrature points per occupied parent
tetrahedron.

If the density basis is already fine enough, `quadrature_subdivision = 1L` is a
reasonable default.

## 2D And 3D Lattices

For `basis = "2d"`, the density basis is a triangular lattice with mesh size
$s_{\text{density}}$.

For `basis = "3d"`, the density basis is a BCC lattice with mesh size
$s_{\text{density}}$.

In both cases, the same origin is used as the parent segmentation basis. New
segmentation fits store this origin in `fit$parameters$origin`. By default,
`spatial_basis_segmentation(origin = NULL)` uses the per-axis median of the
filtered transcript coordinates.

## Regularization

The fitted log-density coefficients are regularized along density-basis graph
edges. For each edge $(a,b)$ with physical length $d_{ab}$:

$$
r_{ab} = \frac{w_a - w_b}{d_{ab}}.
$$

The slope penalty is:

$$
R_{\nabla}(w)
= \lambda
  \frac{1}{|E|}
  \sum_{(a,b)\in E} \psi(r_{ab}).
$$

The available choices for $\psi$ are the same as in
`spatial_basis_log_density_field()`:

- `regularization = "quadratic"`
- `regularization = "huber"`
- `regularization = "bounded"`

An optional graph Laplacian penalty is also available:

$$
R_{\Delta}(w)
= \frac{\lambda_{\Delta}}{2}
  \frac{1}{M}
  \sum_a (\Delta w_a)^2.
$$

Here $\Delta w_a$ is the deviation of $w_a$ from a distance-weighted average of
neighboring density-basis vertices.

The full fitted objective is:

$$
\mathcal{J}(w)
= \frac{1}{n}\mathcal{L}_{\text{data}}
  + R_{\nabla}(w)
  + R_{\Delta}(w).
$$

## Outputs

The returned object contains:

- `total_density`: fitted total density at transcript locations
- `eta`: fitted log-density at transcript locations
- `total_density_basis`: fitted total density at density-basis vertices
- `eta_basis`: fitted log-density coefficients at density-basis vertices
- `basis_points`: coordinates of density-basis vertices
- `basis_lattice`: lattice coordinates of density-basis vertices
- `basis_edges`: density-basis graph used for regularization
- `parent_basis_points`: coordinates of parent segmentation-basis vertices
- `parent_basis_lattice`: lattice coordinates of parent segmentation-basis vertices
- `quadrature`: quadrature basis IDs, weights, and integration weights
- `optim`: the `stats::optim()` result
- `parameters`: effective parameters used for the fit

For basis-level watershed, the most important outputs are:

```r
density_fit$total_density_basis
density_fit$eta_basis
density_fit$basis_points
density_fit$basis_edges
```

These provide density values on a much smaller graph than the transcript graph.

## Example

```r
fit <- spatial_basis_segmentation(
    transcripts_df = transcripts_df,
    cell_signatures = cell_signatures,
    basis = "3d",
    s = 3
)

density_fit <- spatial_basis_total_density_field(
    segmentation_fit = fit,
    basis_subdivision = 2,
    quadrature_subdivision = 1,
    lambda = 1,
    lambda_laplacian = 1,
    regularization = "bounded",
    maxit = 100
)

str(density_fit$total_density_basis)
```

The next step for basis-level partitioning is to evaluate or interpolate the
cell-type spatial prior on the density basis and form:

$$
\rho_k(a) = \rho(a) p_k(a),
$$

where $a$ indexes density-basis vertices. A graph watershed can then be run on
`density_fit$basis_edges` separately for each cell type.
