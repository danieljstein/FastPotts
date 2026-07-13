# Piecewise-linear density field

`spatial_basis_linear_density_field()` fits a total transcript density on the
same triangular or BCC barycentric basis used by `spatial_basis_segmentation()`.
It is parallel to `spatial_basis_log_density_field()`, but the interpolated
quantity is the positive density/intensity itself rather than the log-density.

## Model

For basis point `m`, let

```text
rho_m = exp(alpha_m)
```

and let `phi_m(x)` be the local barycentric interpolation weight at location
`x`. The fitted density is

```text
rho(x) = sum_m phi_m(x) rho_m.
```

At each transcript, only the vertices of its containing simplex have nonzero
weights, so evaluation remains sparse.

The Poisson point-process objective minimized by the implementation is

```text
- sum_i log(sum_m phi_im exp(alpha_m))
+ integral_D sum_m phi_m(x) exp(alpha_m) dx
+ regularization(alpha).
```

This is different from `spatial_basis_log_density_field()`, which uses

```text
rho(x) = exp(sum_m phi_m(x) alpha_m).
```

The log-density model geometrically averages neighboring vertex intensities.
The linear-density model arithmetically averages neighboring vertex
intensities.

## Exact Integration

Because `rho(x)` is piecewise-linear, integration over the occupied simplex
domain is exact. For a simplex `s` with volume `V_s` and `d + 1` vertices,

```text
integral_s rho(x) dx
= V_s / (d + 1) * sum_{m in s} rho_m.
```

Equivalently, if

```text
a_m = integral_D phi_m(x) dx,
```

then

```text
integral_D rho(x) dx = sum_m a_m rho_m.
```

The implementation accumulates this contribution directly from the simplex
volumes in the integration domain. No quadrature points are needed.

## Regularization

Regularization is applied in log-density space:

```text
alpha_m = log(rho_m).
```

The `lambda` penalty smooths neighboring log-density slopes using the same
`quadratic`, `huber`, or `bounded` options as `spatial_basis_log_density_field()`.
The `lambda_laplacian` penalty also acts on `alpha`, penalizing graph-Laplacian
curvature of the log-density field. This makes the penalties act on relative
density changes rather than absolute density differences.

## Output

The returned object follows the compact density-fit layout:

- `total_density_basis`: fitted positive basis-point intensities `rho_m`
- `eta_basis`: fitted log-intensities `alpha_m`
- `domain`: exact occupied-simplex integration domain
- `transcript_basis_id` and `transcript_basis_weight`: sparse interpolation
  design for the fitted transcripts
- `parameters$mode = "linear_density"`
- `parameters$integration = "exact_linear_simplex"`

If `posterior` is supplied, it is not used in the objective. It is only used to
decompose the fitted total density into transcript-by-cell-type densities when
`return_density = TRUE`, matching `spatial_basis_log_density_field()`.
