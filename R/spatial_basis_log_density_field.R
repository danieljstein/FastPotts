#' Fit a log-linear spatial transcript density field
#'
#' Fits a simplified inhomogeneous Poisson point-process density model on the
#' same triangular or BCC barycentric basis used by
#' [spatial_basis_segmentation()]. The model has no boundary attenuation terms:
#'
#' \deqn{
#' \rho(x) = \exp(\eta(x)), \quad
#' \eta(x) = \sum_a \lambda_a(x) w_a.
#' }
#'
#' The integral term is approximated by simplex-local quadrature over occupied
#' simplexes. If `posterior` is supplied, transcript-level cell-type-specific
#' densities are returned as `total_density * posterior`.
#'
#' @param transcripts_df Transcript-level data frame.
#' @param posterior Optional numeric transcript-by-cell-type posterior matrix.
#'   If supplied, rows are normalized and used only to decompose the fitted
#'   total density into cell-type-specific densities.
#' @param return_density Logical; if `TRUE` and `posterior` is supplied, return
#'   the transcript-by-cell-type density matrix `total_density * posterior`.
#'   Set to `FALSE` for large datasets when only `total_density` is needed.
#' @param basis Character; `"2d"`/`"tri"` or `"3d"`/`"bcc"`.
#' @param s Positive spatial basis mesh size.
#' @param x,y,z Character coordinate column names.
#' @param origin Optional lattice origin.
#' @param quadrature_subdivision Positive integer subdivision factor for
#'   occupied-simplex quadrature.
#' @param store_quadrature_coords Logical; if `TRUE`, store quadrature point
#'   coordinates in the returned object. The fitting objective does not need
#'   these coordinates, so the default `FALSE` is more memory efficient for
#'   large 3D datasets.
#' @param lambda Non-negative smoothing strength on neighboring basis-point
#'   log-density slopes.
#' @param regularization Character; one of `"quadratic"`, `"huber"`, or
#'   `"bounded"`.
#' @param delta Positive Huber transition scale. Used only when
#'   `regularization = "huber"`.
#' @param sigma Positive bounded-penalty slope scale. Used only when
#'   `regularization = "bounded"`.
#' @param lambda_laplacian Non-negative strength for a graph Laplacian
#'   curvature penalty on the basis-point log-density field. This penalizes
#'   deviations from a distance-weighted neighbor average.
#' @param maxit Maximum L-BFGS iterations.
#' @param reltol Relative convergence tolerance.
#' @param n_threads Integer number of OpenMP threads. If `NULL`, uses runtime
#'   default.
#' @param show_progress Logical; print progress messages.
#'
#' @return A list with fitted total density, optional cell-type density,
#'   fitted basis coefficients, quadrature design, basis design, optimizer
#'   result, and effective parameters.
#' @export
spatial_basis_log_density_field <- function(
    transcripts_df,
    posterior = NULL,
    return_density = !is.null(posterior),
    basis = c("3d", "2d", "tri", "bcc"),
    s,
    x = "x_location",
    y = "y_location",
    z = "z_location",
    origin = NULL,
    quadrature_subdivision = 4L,
    store_quadrature_coords = FALSE,
    lambda = 0.1,
    regularization = c("bounded", "quadratic", "huber"),
    delta = 1,
    sigma = 1,
    lambda_laplacian = 0,
    maxit = 100L,
    reltol = 1e-6,
    n_threads = NULL,
    show_progress = TRUE
) {
    basis = normalize_spatial_basis(basis)
    regularization = match.arg(regularization)
    lattice_basis = if (basis == "2d") "tri" else "bcc"
    d = if (basis == "2d") 2L else 3L
    coord_cols = if (basis == "2d") c(x, y) else c(x, y, z)
    missing_cols = setdiff(coord_cols, colnames(transcripts_df))
    if (length(missing_cols) > 0L) {
        stop("Missing coordinate column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    if (is.null(origin)) {
        origin = rep(0, d)
    }
    if (length(origin) != d || any(!is.finite(origin))) {
        stop("origin must be a finite numeric vector with length matching the selected basis.", call. = FALSE)
    }
    if (length(s) != 1L || !is.finite(s) || s <= 0) {
        stop("s must be a positive finite scalar.", call. = FALSE)
    }
    if (length(quadrature_subdivision) != 1L ||
        !is.finite(quadrature_subdivision) ||
        quadrature_subdivision < 1L) {
        stop("quadrature_subdivision must be a positive integer.", call. = FALSE)
    }
    quadrature_subdivision = as.integer(quadrature_subdivision)
    if (!is.logical(store_quadrature_coords) || length(store_quadrature_coords) != 1L || is.na(store_quadrature_coords)) {
        stop("store_quadrature_coords must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(return_density) || length(return_density) != 1L || is.na(return_density)) {
        stop("return_density must be TRUE or FALSE.", call. = FALSE)
    }
    if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
        stop("lambda must be a non-negative finite scalar.", call. = FALSE)
    }
    if (length(delta) != 1L || !is.finite(delta) || delta <= 0) {
        stop("delta must be a positive finite scalar.", call. = FALSE)
    }
    if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
        stop("sigma must be a positive finite scalar.", call. = FALSE)
    }
    if (length(lambda_laplacian) != 1L || !is.finite(lambda_laplacian) || lambda_laplacian < 0) {
        stop("lambda_laplacian must be a non-negative finite scalar.", call. = FALSE)
    }
    if (maxit < 1L) {
        stop("maxit must be positive.", call. = FALSE)
    }
    if (is.null(n_threads)) {
        n_threads = 0L
    } else if (length(n_threads) != 1L || !is.finite(n_threads) || n_threads < 1) {
        stop("n_threads must be NULL or a positive integer.", call. = FALSE)
    } else {
        n_threads = as.integer(n_threads)
    }
    basis_n_threads = if (n_threads == 0L) NULL else n_threads

    coords = as.matrix(transcripts_df[, coord_cols, drop = FALSE])
    storage.mode(coords) = "double"
    if (any(!is.finite(coords))) {
        stop("Coordinate columns must contain finite values.", call. = FALSE)
    }

    if (!is.null(posterior)) {
        posterior = normalize_posterior_matrix(posterior, nrow(transcripts_df))
        cell_types = colnames(posterior)
        if (is.null(cell_types)) {
            cell_types = paste0("type", seq_len(ncol(posterior)))
            colnames(posterior) = cell_types
        }
    } else {
        cell_types = character()
    }

    if (show_progress) {
        message("Computing spatial basis interpolation...")
    }
    bary = if (basis == "2d") {
        tri_barycentric(coords, s = s, origin = origin, n_threads = basis_n_threads)
    } else {
        bcc_barycentric(coords, s = s, origin = origin, n_threads = basis_n_threads)
    }
    design = basis_design_from_barycentric(bary)

    if (show_progress) {
        message("Building occupied-simplex quadrature...")
    }
    quadrature = make_density_quadrature_simplex(
        design = design,
        d = d,
        subdivision = quadrature_subdivision,
        store_coords = store_quadrature_coords
    )

    if (show_progress) {
        message("Building basis edge design...")
    }
    basis_edges = build_lattice_neighbor_edges(design$basis_lattice, basis = lattice_basis, s = s)

    obs_basis_id = design$basis_id - 1L
    obs_basis_weight = design$basis_weight
    quad_basis_id = quadrature$basis_id - 1L
    quad_basis_weight = quadrature$basis_weight

    volume = sum(quadrature$weight)
    par0 = rep(log(pmax(nrow(coords) / volume, 1e-8)), nrow(design$basis_lattice))
    regularization_id = match(regularization, c("quadratic", "huber", "bounded")) - 1L

    objective = function(par) {
        spatial_log_density_objective_cpp(
            par = par,
            obs_basis_id = obs_basis_id,
            obs_basis_weight = obs_basis_weight,
            quad_basis_id = quad_basis_id,
            quad_basis_weight = quad_basis_weight,
            quad_weight = quadrature$weight,
            edge_from = as.integer(basis_edges[, "from"] - 1L),
            edge_to = as.integer(basis_edges[, "to"] - 1L),
            edge_distance = as.numeric(basis_edges[, "distance"]),
            lambda = lambda,
            regularization = regularization_id,
            delta = delta,
            sigma = sigma,
            lambda_laplacian = lambda_laplacian,
            n_threads = n_threads
        )
    }

    if (show_progress) {
        message("Optimizing log-density field...")
    }
    opt = stats::optim(
        par = par0,
        fn = function(par) objective(par)$value,
        gr = function(par) objective(par)$gradient,
        method = "L-BFGS-B",
        control = list(maxit = as.integer(maxit), factr = reltol / .Machine$double.eps)
    )
    warn_spatial_basis_optim_status(opt, maxit)

    pred = spatial_log_density_predict_cpp(
        par = opt$par,
        basis_id = obs_basis_id,
        basis_weight = obs_basis_weight,
        n_threads = n_threads
    )
    total_density = as.numeric(pred$density)
    eta = as.numeric(pred$eta)

    if (!is.null(posterior) && return_density) {
        density = posterior * total_density
        colnames(density) = cell_types
    } else {
        density = NULL
    }

    list(
        density = density,
        total_density = total_density,
        eta = eta,
        eta_basis = opt$par,
        posterior = posterior,
        basis_points = design$basis_points,
        basis_lattice = design$basis_lattice,
        basis_edges = basis_edges,
        quadrature = quadrature,
        optim = opt,
        parameters = list(
            basis = basis,
            s = s,
            origin = origin,
            quadrature_subdivision = quadrature_subdivision,
            store_quadrature_coords = store_quadrature_coords,
            return_density = return_density,
            lambda = lambda,
            regularization = regularization,
            delta = delta,
            sigma = sigma,
            lambda_laplacian = lambda_laplacian,
            maxit = maxit,
            reltol = reltol
        )
    )
}
