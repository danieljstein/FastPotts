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

new_basis_design_accumulator <- function(d) {
    list(
        key_to_id = integer(),
        basis_lattice = matrix(integer(), nrow = 0L, ncol = d),
        basis_points = matrix(numeric(), nrow = 0L, ncol = d)
    )
}

append_basis_design_chunk <- function(accumulator, design) {
    local_key = do.call(paste, c(as.data.frame(design$basis_lattice), sep = ":"))
    global_id = unname(accumulator$key_to_id[local_key])
    missing = is.na(global_id)
    if (any(missing)) {
        new_id = seq.int(
            length(accumulator$key_to_id) + 1L,
            length(accumulator$key_to_id) + sum(missing)
        )
        names(new_id) = local_key[missing]
        accumulator$key_to_id = c(accumulator$key_to_id, new_id)
        accumulator$basis_lattice = rbind(
            accumulator$basis_lattice,
            design$basis_lattice[missing, , drop = FALSE]
        )
        accumulator$basis_points = rbind(
            accumulator$basis_points,
            design$basis_points[missing, , drop = FALSE]
        )
        global_id[missing] = new_id
    }

    design$basis_id[] = global_id[design$basis_id]
    list(
        accumulator = accumulator,
        basis_id = design$basis_id,
        basis_weight = design$basis_weight
    )
}

compute_density_basis_design_chunked <- function(
    coords,
    basis,
    s,
    origin,
    n_threads,
    chunk_size,
    accumulator = NULL
) {
    coords = as.matrix(coords)
    d = ncol(coords)
    n = nrow(coords)
    if (is.null(accumulator)) {
        accumulator = new_basis_design_accumulator(d)
    }
    n_active = if (basis == "2d") 3L else 4L
    basis_id = matrix(NA_integer_, nrow = n, ncol = n_active)
    basis_weight = matrix(NA_real_, nrow = n, ncol = n_active)
    if (n == 0L) {
        return(list(
            accumulator = accumulator,
            basis_id = basis_id,
            basis_weight = basis_weight
        ))
    }

    chunk_starts = seq.int(1L, n, by = chunk_size)
    for (start in chunk_starts) {
        end = min(start + chunk_size - 1L, n)
        chunk_coords = coords[start:end, , drop = FALSE]
        bary = if (basis == "2d") {
            tri_barycentric(chunk_coords, s = s, origin = origin, n_threads = n_threads)
        } else {
            bcc_barycentric(chunk_coords, s = s, origin = origin, n_threads = n_threads)
        }
        chunk_design = basis_design_from_barycentric(bary)
        merged = append_basis_design_chunk(accumulator, chunk_design)
        accumulator = merged$accumulator
        basis_id[start:end, ] = merged$basis_id
        basis_weight[start:end, ] = merged$basis_weight
    }

    list(
        accumulator = accumulator,
        basis_id = basis_id,
        basis_weight = basis_weight
    )
}

align_density_basis_design_to_parent <- function(obs_basis_id, quad_basis_id, basis_lattice_density, basis_lattice_parent) {
    parent_key = do.call(paste, c(as.data.frame(basis_lattice_parent), sep = ":"))
    density_key = do.call(paste, c(as.data.frame(basis_lattice_density), sep = ":"))
    parent_lookup = seq_along(parent_key)
    names(parent_lookup) = parent_key

    matched = unname(parent_lookup[density_key])
    missing = is.na(matched)
    if (any(missing)) {
        stop(
            "The density basis does not align with segmentation_fit$basis_lattice; ",
            sum(missing),
            " active basis point(s) were not found. ",
            "Check that basis, s, origin, and transcripts_df match the segmentation fit.",
            call. = FALSE
        )
    }

    list(
        obs_basis_id = matrix(matched[obs_basis_id], nrow = nrow(obs_basis_id), ncol = ncol(obs_basis_id)),
        quad_basis_id = matrix(matched[quad_basis_id], nrow = nrow(quad_basis_id), ncol = ncol(quad_basis_id))
    )
}

finalize_density_basis_design <- function(obs_basis_id, quad_basis_id, basis_lattice, basis_points) {
    key = do.call(paste, c(as.data.frame(basis_lattice), sep = ":"))
    order_id = order(key)
    old_to_new = integer(length(order_id))
    old_to_new[order_id] = seq_along(order_id)

    list(
        obs_basis_id = matrix(old_to_new[obs_basis_id], nrow = nrow(obs_basis_id), ncol = ncol(obs_basis_id)),
        quad_basis_id = matrix(old_to_new[quad_basis_id], nrow = nrow(quad_basis_id), ncol = ncol(quad_basis_id)),
        basis_lattice = basis_lattice[order_id, , drop = FALSE],
        basis_points = basis_points[order_id, , drop = FALSE]
    )
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
#' @param interpolation_chunk_size Positive integer number of transcript or
#'   quadrature coordinates to pass to the barycentric interpolation backend at
#'   once. Smaller values reduce peak memory during density-basis interpolation.
#' @param lambda Non-negative smoothing strength on neighboring basis-point
#'   log-density slopes.
#' @param regularization Character; one of `"quadratic"`, `"huber"`, or
#'   `"bounded"`.
#' @param delta Positive Huber transition scale. Used only when
#'   `regularization = "huber"`.
#' @param sigma Positive bounded-penalty slope scale. Used only when
#'   `regularization = "bounded"`.
#' @param lambda_laplacian Non-negative strength for a graph Laplacian
#'   curvature penalty on the basis-point log-density field. The supplied value
#'   is scaled internally by `density_s^4`, where `density_s = s /
#'   basis_subdivision`, so that the parameter is approximately invariant to the
#'   density mesh size.
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
    interpolation_chunk_size = 250000L,
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
    if (length(interpolation_chunk_size) != 1L ||
        !is.finite(interpolation_chunk_size) ||
        interpolation_chunk_size < 1L) {
        stop("interpolation_chunk_size must be a positive integer.", call. = FALSE)
    }
    interpolation_chunk_size = as.integer(interpolation_chunk_size)
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
    obs_density_design = compute_density_basis_design_chunked(
        coords = coords,
        basis = basis,
        s = density_s,
        origin = origin,
        n_threads = basis_n_threads,
        chunk_size = interpolation_chunk_size
    )
    quad_density_design = compute_density_basis_design_chunked(
        coords = parent_quadrature$coords,
        basis = basis,
        s = density_s,
        origin = origin,
        n_threads = basis_n_threads,
        chunk_size = interpolation_chunk_size,
        accumulator = obs_density_design$accumulator
    )

    finalized_density_design = finalize_density_basis_design(
        obs_basis_id = obs_density_design$basis_id,
        quad_basis_id = quad_density_design$basis_id,
        basis_lattice = quad_density_design$accumulator$basis_lattice,
        basis_points = quad_density_design$accumulator$basis_points
    )
    obs_density_design$basis_id = finalized_density_design$obs_basis_id
    quad_density_design$basis_id = finalized_density_design$quad_basis_id
    basis_lattice_density = finalized_density_design$basis_lattice
    basis_points_density = finalized_density_design$basis_points

    if (basis_subdivision == 1L) {
        aligned_density_design = align_density_basis_design_to_parent(
            obs_basis_id = obs_density_design$basis_id,
            quad_basis_id = quad_density_design$basis_id,
            basis_lattice_density = basis_lattice_density,
            basis_lattice_parent = basis_lattice
        )
        obs_density_design$basis_id = aligned_density_design$obs_basis_id
        quad_density_design$basis_id = aligned_density_design$quad_basis_id
        basis_lattice_density = basis_lattice
        basis_points_density = basis_points
    }
    basis_edges = if (basis_subdivision == 1L) segmentation_fit$basis_edges else NULL
    if (is.null(basis_edges)) {
        if (show_progress) {
            message("Building density basis edge design...")
        }
        basis_edges = build_lattice_neighbor_edges(basis_lattice_density, basis = lattice_basis, s = density_s)
    }
    lambda_laplacian_effective = lambda_laplacian / (density_s^4)
    edge_from = as.integer(basis_edges[, "from"] - 1L)
    edge_to = as.integer(basis_edges[, "to"] - 1L)
    edge_distance = as.numeric(basis_edges[, "distance"])

    quadrature = list(
        coords = if (isTRUE(store_quadrature_coords)) parent_quadrature$coords else NULL,
        weight = parent_quadrature$weight,
        basis_id = quad_density_design$basis_id,
        basis_weight = quad_density_design$basis_weight,
        bounds = if (isTRUE(store_quadrature_coords)) parent_quadrature$bounds else NULL,
        method = "parent_simplex_density_basis",
        subdivision = as.integer(quadrature_subdivision),
        basis_subdivision = as.integer(basis_subdivision),
        effective_subdivision = effective_quadrature_subdivision,
        n_simplex = parent_quadrature$n_simplex
    )

    obs_basis_id = obs_density_design$basis_id - 1L
    obs_basis_weight = obs_density_design$basis_weight
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
            edge_from = edge_from,
            edge_to = edge_to,
            edge_distance = edge_distance,
            lambda = lambda,
            regularization = regularization_id,
            delta = delta,
            sigma = sigma,
            lambda_laplacian = lambda_laplacian_effective,
            n_threads = n_threads
        )
    }

    if (show_progress) {
        message("Optimizing total log-density field...")
    }
    cached = make_cached_spatial_basis_objective(objective)
    opt = stats::optim(
        par = par0,
        fn = cached$fn,
        gr = cached$gr,
        method = "L-BFGS-B",
        control = list(maxit = as.integer(maxit), factr = reltol / .Machine$double.eps)
    )
    attr(opt, "objective_evaluations") = cached$n_eval()
    attr(opt, "objective_cache_hits") = cached$n_hit()
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
        transcript_basis_id = obs_density_design$basis_id,
        transcript_basis_weight = obs_density_design$basis_weight,
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
            interpolation_chunk_size = interpolation_chunk_size,
            lambda = lambda,
            regularization = regularization,
            delta = delta,
            sigma = sigma,
            lambda_laplacian = lambda_laplacian,
            lambda_laplacian_effective = lambda_laplacian_effective,
            maxit = maxit,
            reltol = reltol
        )
    )
}
