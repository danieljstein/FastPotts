#' Barycentric coordinates on a body-centered cubic lattice
#'
#' Computes the four active piecewise-linear basis functions for each point on
#' a body-centered cubic (BCC) lattice. The lattice is the union of
#' `(2*s*i, 2*s*j, 2*s*k)` and `(s + 2*s*i, s + 2*s*j, s + 2*s*k)`, optionally
#' shifted by `origin`.
#'
#' For each query point, the function finds the containing BCC Delaunay
#' tetrahedron and returns its four vertices and barycentric weights. Away from
#' tetrahedron boundaries exactly four basis functions are non-zero. On faces,
#' edges, or vertices, some returned weights may be zero.
#'
#' @param coords Numeric matrix or data frame with three columns containing
#'   query coordinates.
#' @param s Positive numeric mesh size.
#' @param origin Numeric vector of length three giving the lattice origin.
#' @param tol Numeric tolerance used for tetrahedron boundary checks.
#'
#' @return A list with components:
#' \describe{
#'   \item{weights}{Numeric matrix with four barycentric weights per row of
#'     `coords`.}
#'   \item{points}{Numeric array of dimension `n x 4 x 3` containing the four
#'     active lattice points for each query point.}
#'   \item{lattice}{Integer array of dimension `n x 4 x 3` containing the
#'     normalized BCC lattice coordinates before multiplying by `s` and adding
#'     `origin`.}
#' }
#'
#' @examples
#' coords <- matrix(c(0.75, 0.75, 0.25), ncol = 3)
#' bcc <- bcc_barycentric(coords, s = 1)
#' rowSums(bcc$weights)
#'
#' @export
bcc_barycentric <- function(coords, s, origin = c(0, 0, 0), tol = 1e-10) {
    coords <- as.matrix(coords)

    if (ncol(coords) != 3L) {
        stop("coords must have exactly 3 columns.")
    }

    if (
        length(s) != 1L ||
        !is.finite(s) ||
        s <= 0
    ) {
        stop("s must be a positive finite number.")
    }

    if (length(origin) != 3L || any(!is.finite(origin))) {
        stop("origin must be a finite numeric vector of length 3.")
    }

    if (
        length(tol) != 1L ||
        !is.finite(tol) ||
        tol < 0
    ) {
        stop("tol must be a non-negative finite number.")
    }

    bcc_barycentric_cpp(coords, s = as.numeric(s), origin = as.numeric(origin), tol = as.numeric(tol))
}
