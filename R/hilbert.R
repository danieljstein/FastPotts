#' Compute Hilbert indices for 2D or 3D coordinates, scaled to the unit cube.
#' 
#' @param coords A numeric matrix with 2 or 3 columns, containing the coordinates
#'   for which to compute Hilbert indices.
#' @param level The number of bits to use for each coordinate. If `NULL`, uses the maximum possible level that
#'   can be represented as a double (i.e. `floor(53 / d)` where `d` is the number of dimensions).
#' @param return_scaled If `TRUE`, returns a list with components `hilbert_index` and `coords_scaled`,
#'   where `hilbert_index` is the Hilbert index computed from the scaled coordinates,
#'   and `coords_scaled` is the input coordinates scaled to the unit cube according to the specified level.
#'   If `FALSE`, returns just the Hilbert index.
#' 
#' @export
hilbert_index <- function(coords, level = NULL, return_scaled = FALSE) {
    coords <- as.matrix(coords)
    d <- ncol(coords)

    if (!d %in% c(2L, 3L)) {
        stop("coords must have 2 or 3 columns.")
    }

    if (is.null(level)) {
        level <- floor(53 / d)
    }

    if (
        length(level) != 1L ||
        !is.finite(level) ||
        level < 1 ||
        level != as.integer(level)
    ) {
        stop("level must be a positive integer.")
    }

    level <- as.integer(level)

    if (level > 30L) {
        stop("level must be <= 30 because scaled coordinates are stored as R integers.")
    }

    if (d * level > 53) {
        warning("d * level > 53, so returned Hilbert indices may not be exactly representable as doubles.")
    }

    hilbert_index_scaled_cpp(coords, bits = level, return_scaled = return_scaled)
}
