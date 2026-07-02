infer_spatial_basis_mesh_size <- function(segmentation_fit, basis) {
    if (is.null(segmentation_fit$basis_edges) || !"distance" %in% colnames(segmentation_fit$basis_edges)) {
        stop("s must be supplied when segmentation_fit$basis_edges is unavailable.", call. = FALSE)
    }
    distance = as.numeric(segmentation_fit$basis_edges[, "distance"])
    distance = distance[is.finite(distance) & distance > 0]
    if (length(distance) == 0L) {
        stop("Could not infer s from segmentation_fit$basis_edges.", call. = FALSE)
    }
    if (basis == "2d") {
        return(stats::median(distance))
    }
    max(distance)
}

infer_spatial_basis_origin <- function(basis_lattice, basis_points, basis, s) {
    if (nrow(basis_lattice) < 1L) {
        stop("basis_lattice must contain at least one row.", call. = FALSE)
    }
    if (basis == "2d") {
        ij = basis_lattice[1L, ]
        return(c(
            basis_points[1L, 1L] - s * (ij[1L] + 0.5 * ij[2L]),
            basis_points[1L, 2L] - s * sqrt(3) * ij[2L] / 2
        ))
    }
    basis_points[1L, ] - (s / 2) * basis_lattice[1L, ]
}

align_design_to_basis <- function(design, basis_lattice, basis_points) {
    basis_key = do.call(paste, c(as.data.frame(basis_lattice), sep = ":"))
    design_key = do.call(paste, c(as.data.frame(design$basis_lattice), sep = ":"))
    basis_lookup = seq_along(basis_key)
    names(basis_lookup) = basis_key

    matched = unname(basis_lookup[design_key])
    missing = is.na(matched)
    if (any(missing)) {
        stop(
            "The transcript coordinates do not align with the supplied basis; ",
            sum(missing),
            " active basis point(s) were not found in segmentation_fit$basis_lattice. ",
            "Check that basis, s, origin, and transcripts_df match the segmentation fit.",
            call. = FALSE
        )
    }

    remap = seq_along(matched)
    names(remap) = remap
    remap[] = matched
    design$basis_id[] = unname(remap[as.character(design$basis_id)])
    design$basis_lattice = basis_lattice
    design$basis_points = basis_points
    design
}

#' Fit total transcript density on an existing spatial basis
#'
#' Fits the same log-linear inhomogeneous Poisson density model as
#' [spatial_basis_log_density_field()], but reuses the basis points from an
#' existing [spatial_basis_segmentation()] result. The model ignores cell types
#' and estimates one total transcript density field:
#'
#' \deqn{
#' \rho(x) = \exp(\eta(x)), \quad
#' \eta(x) = \sum_a \lambda_a(x) w_a.
#' }
#'
#' @param segmentation_fit Result from [spatial_basis_segmentation()] containing
#'   `basis_lattice`, `basis_points`, and preferably `basis_edges`.
#' @param transcripts_df Transcript-level data frame. Defaults to
#'   `segmentation_fit$transcripts_df`.
#' @param basis Character; `"2d"`/`"tri"` or `"3d"`/`"bcc"`. If `NULL`, inferred
#'   from the number of basis-point coordinate columns.
#' @param s Positive spatial basis mesh size. If `NULL`, inferred from
#'   `segmentation_fit$basis_edges`.
#' @param x,y,z Character coordinate column names.
#' @param origin Optional lattice origin. If `NULL`, inferred from
#'   `basis_lattice`, `basis_points`, and `s`.
#' @param basis_subdivision Positive integer refinement factor for the density
#'   basis relative to the segmentation basis. The density mesh size is
#'   `s / basis_subdivision`. The integration domain is still defined by the
#'   occupied simplexes of the segmentation basis, so child regions with no
#'   transcripts still contribute to the Poisson integral.
#' @param quadrature_subdivision Positive integer subdivision factor for
#'   quadrature within each density-basis subdivision. The effective quadrature
#'   subdivision of each occupied segmentation simplex is
#'   `basis_subdivision * quadrature_subdivision`.
#' @param store_quadrature_coords Logical; if `TRUE`, store quadrature point
#'   coordinates in the returned object. The fitting objective does not need
#'   these coordinates, so the default `FALSE` is more memory efficient.
#' @param lambda Non-negative smoothing strength on neighboring basis-point
#'   log-density slopes.
#' @param regularization Character; one of `"quadratic"`, `"huber"`, or
#'   `"bounded"`.
#' @param delta Positive Huber transition scale. Used only when
#'   `regularization = "huber"`.
#' @param sigma Positive bounded-penalty slope scale. Used only when
#'   `regularization = "bounded"`.
#' @param lambda_laplacian Non-negative strength for a graph Laplacian
#'   curvature penalty on the basis-point log-density field.
#' @param maxit Maximum L-BFGS iterations.
#' @param reltol Relative convergence tolerance.
#' @param n_threads Integer number of OpenMP threads. If `NULL`, uses runtime
#'   default.
#' @param show_progress Logical; print progress messages.
#'
#' @return A list with fitted transcript-level and basis-point total densities,
#'   fitted basis coefficients, quadrature design, reused basis design,
#'   optimizer result, and effective parameters.
#' @export
spatial_basis_total_density_field <- function(
    segmentation_fit,
    transcripts_df = segmentation_fit$transcripts_df,
    basis = NULL,
    s = NULL,
    x = "x_location",
    y = "y_location",
    z = "z_location",
    origin = NULL,
    basis_subdivision = 1L,
    quadrature_subdivision = 1L,
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
    if (is.null(segmentation_fit$basis_lattice) || is.null(segmentation_fit$basis_points)) {
        stop("segmentation_fit must contain basis_lattice and basis_points.", call. = FALSE)
    }
    basis_lattice = as.matrix(segmentation_fit$basis_lattice)
    basis_points = as.matrix(segmentation_fit$basis_points)
    storage.mode(basis_lattice) = "integer"
    storage.mode(basis_points) = "double"
    if (nrow(basis_lattice) != nrow(basis_points)) {
        stop("segmentation_fit$basis_lattice and basis_points must have the same number of rows.", call. = FALSE)
    }
    if (is.null(basis)) {
        basis = if (!is.null(segmentation_fit$parameters$basis)) {
            segmentation_fit$parameters$basis
        } else if (ncol(basis_points) == 2L) {
            "2d"
        } else if (ncol(basis_points) == 3L) {
            "3d"
        } else {
            NA_character_
        }
    }
    basis = normalize_spatial_basis(basis)
    regularization = match.arg(regularization)
    lattice_basis = if (basis == "2d") "tri" else "bcc"
    d = if (basis == "2d") 2L else 3L
    if (ncol(basis_lattice) != d || ncol(basis_points) != d) {
        stop("Supplied basis geometry does not match selected basis dimensionality.", call. = FALSE)
    }
    if (is.null(s)) {
        s = if (!is.null(segmentation_fit$parameters$s)) {
            segmentation_fit$parameters$s
        } else {
            infer_spatial_basis_mesh_size(segmentation_fit, basis)
        }
    }
    if (length(s) != 1L || !is.finite(s) || s <= 0) {
        stop("s must be a positive finite scalar.", call. = FALSE)
    }
    if (is.null(origin)) {
        origin = if (!is.null(segmentation_fit$parameters$origin)) {
            segmentation_fit$parameters$origin
        } else {
            infer_spatial_basis_origin(basis_lattice, basis_points, basis, s)
        }
    }
    if (length(origin) != d || any(!is.finite(origin))) {
        stop("origin must be a finite numeric vector with length matching the selected basis.", call. = FALSE)
    }
    if (length(quadrature_subdivision) != 1L ||
        !is.finite(quadrature_subdivision) ||
        quadrature_subdivision < 1L) {
        stop("quadrature_subdivision must be a positive integer.", call. = FALSE)
    }
    if (length(basis_subdivision) != 1L ||
        !is.finite(basis_subdivision) ||
        basis_subdivision < 1L) {
        stop("basis_subdivision must be a positive integer.", call. = FALSE)
    }
    basis_subdivision = as.integer(basis_subdivision)
    quadrature_subdivision = as.integer(quadrature_subdivision)
    if (!is.logical(store_quadrature_coords) || length(store_quadrature_coords) != 1L || is.na(store_quadrature_coords)) {
        stop("store_quadrature_coords must be TRUE or FALSE.", call. = FALSE)
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

    coord_cols = if (basis == "2d") c(x, y) else c(x, y, z)
    missing_cols = setdiff(coord_cols, colnames(transcripts_df))
    if (length(missing_cols) > 0L) {
        stop("Missing coordinate column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    coords = as.matrix(transcripts_df[, coord_cols, drop = FALSE])
    storage.mode(coords) = "double"
    if (any(!is.finite(coords))) {
        stop("Coordinate columns must contain finite values.", call. = FALSE)
    }

    if (show_progress) {
        message("Computing transcript interpolation on parent basis...")
    }
    parent_bary = if (basis == "2d") {
        tri_barycentric(coords, s = s, origin = origin, n_threads = basis_n_threads)
    } else {
        bcc_barycentric(coords, s = s, origin = origin, n_threads = basis_n_threads)
    }
    parent_design = align_design_to_basis(basis_design_from_barycentric(parent_bary), basis_lattice, basis_points)

    if (show_progress) {
        message("Building parent-domain quadrature...")
    }
    effective_quadrature_subdivision = as.integer(basis_subdivision * quadrature_subdivision)
    parent_quadrature = make_density_quadrature_simplex(
        design = parent_design,
        d = d,
        subdivision = effective_quadrature_subdivision,
        store_coords = TRUE
    )

    density_s = s / basis_subdivision
    if (show_progress) {
        message("Computing interpolation on density basis...")
    }
    all_coords = rbind(coords, parent_quadrature$coords)
    density_bary = if (basis == "2d") {
        tri_barycentric(all_coords, s = density_s, origin = origin, n_threads = basis_n_threads)
    } else {
        bcc_barycentric(all_coords, s = density_s, origin = origin, n_threads = basis_n_threads)
    }
    design = basis_design_from_barycentric(density_bary)
    obs_rows = seq_len(nrow(coords))
    quad_rows = seq.int(nrow(coords) + 1L, nrow(all_coords))

    if (basis_subdivision == 1L) {
        design = align_design_to_basis(design, basis_lattice, basis_points)
    }

    basis_lattice_density = design$basis_lattice
    basis_points_density = design$basis_points
    basis_edges = if (basis_subdivision == 1L) segmentation_fit$basis_edges else NULL
    if (is.null(basis_edges)) {
        if (show_progress) {
            message("Building density basis edge design...")
        }
        basis_edges = build_lattice_neighbor_edges(basis_lattice_density, basis = lattice_basis, s = density_s)
    }

    quadrature = list(
        coords = if (isTRUE(store_quadrature_coords)) parent_quadrature$coords else NULL,
        weight = parent_quadrature$weight,
        basis_id = design$basis_id[quad_rows, , drop = FALSE],
        basis_weight = design$basis_weight[quad_rows, , drop = FALSE],
        bounds = if (isTRUE(store_quadrature_coords)) parent_quadrature$bounds else NULL,
        method = "parent_simplex_density_basis",
        subdivision = as.integer(quadrature_subdivision),
        basis_subdivision = as.integer(basis_subdivision),
        effective_subdivision = effective_quadrature_subdivision,
        n_simplex = parent_quadrature$n_simplex
    )

    obs_basis_id = design$basis_id[obs_rows, , drop = FALSE] - 1L
    obs_basis_weight = design$basis_weight[obs_rows, , drop = FALSE]
    quad_basis_id = quadrature$basis_id - 1L
    quad_basis_weight = quadrature$basis_weight

    volume = sum(quadrature$weight)
    par0 = rep(log(pmax(nrow(coords) / volume, 1e-8)), nrow(basis_lattice_density))
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
        message("Optimizing total log-density field...")
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

    list(
        total_density = as.numeric(pred$density),
        eta = as.numeric(pred$eta),
        total_density_basis = as.numeric(exp(opt$par)),
        eta_basis = opt$par,
        basis_points = basis_points_density,
        basis_lattice = basis_lattice_density,
        basis_edges = basis_edges,
        parent_basis_points = basis_points,
        parent_basis_lattice = basis_lattice,
        quadrature = quadrature,
        optim = opt,
        parameters = list(
            basis = basis,
            s = s,
            density_s = density_s,
            origin = origin,
            basis_subdivision = basis_subdivision,
            quadrature_subdivision = quadrature_subdivision,
            effective_quadrature_subdivision = effective_quadrature_subdivision,
            store_quadrature_coords = store_quadrature_coords,
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
