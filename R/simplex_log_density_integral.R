#' Integrate a log-linear density over triangles or tetrahedra
#'
#' Computes the integral and vertex-value gradient for a log-linear density on
#' each simplex:
#'
#' \deqn{
#' I(h) = \int_S \exp\left(\sum_i \lambda_i(x) h_i\right) dx.
#' }
#'
#' The returned gradient column `i` is
#'
#' \deqn{
#' \frac{\partial I}{\partial h_i}
#' =
#' \int_S \lambda_i(x)
#' \exp\left(\sum_j \lambda_j(x) h_j\right) dx.
#' }
#'
#' The implementation uses a Taylor expansion when vertex values are close, a
#' shifted divided-difference formula when they are well separated, and a
#' deterministic Gauss fallback for nearly tied subsets with large spread.
#'
#' @param h Numeric matrix with one row per simplex and either three columns
#'   for triangles or four columns for tetrahedra. Entries are vertex
#'   log-density values.
#' @param volume Numeric scalar or vector with one simplex volume per row of
#'   `h`.
#' @param taylor_radius Non-negative spread threshold for using the Taylor
#'   expansion.
#' @param close_tol Relative pairwise-gap threshold for the Gauss fallback when
#'   the full simplex spread is not small.
#' @param taylor_tol Positive Taylor stopping tolerance.
#' @param taylor_max_terms Maximum Taylor order.
#' @param gauss_order Number of Gauss-Legendre nodes per transformed simplex
#'   dimension for the rare fallback path.
#' @param n_threads Integer number of OpenMP threads. If `NULL`, uses runtime
#'   default.
#'
#' @return A list with `integral`, `gradient`, and integer `method`, where
#'   method is 0 for divided differences, 1 for Taylor, and 2 for Gauss
#'   fallback.
#' @export
simplex_log_density_integral <- function(
    h,
    volume = 1,
    taylor_radius = 1,
    close_tol = 1e-6,
    taylor_tol = 1e-12,
    taylor_max_terms = 80L,
    gauss_order = 20L,
    n_threads = NULL
) {
    h = as.matrix(h)
    storage.mode(h) = "double"
    volume = as.numeric(volume)
    if (is.null(n_threads)) {
        n_threads = 0L
    } else if (length(n_threads) != 1L || !is.finite(n_threads) || n_threads < 1) {
        stop("n_threads must be NULL or a positive integer.", call. = FALSE)
    } else {
        n_threads = as.integer(n_threads)
    }
    simplex_log_density_integral_cpp(
        h = h,
        volume = volume,
        taylor_radius = as.numeric(taylor_radius),
        close_tol = as.numeric(close_tol),
        taylor_tol = as.numeric(taylor_tol),
        taylor_max_terms = as.integer(taylor_max_terms),
        gauss_order = as.integer(gauss_order),
        n_threads = n_threads
    )
}
