softmax_rows <- function(x) {
    x = as.matrix(x)
    row_max = do.call(pmax, as.data.frame(x))
    z = exp(sweep(x, 1L, row_max, "-"))
    sweep(z, 1L, rowSums(z), "/")
}

evaluate_spatial_prior_on_points <- function(
    segmentation_fit,
    coords,
    basis = NULL,
    s = NULL,
    origin = NULL,
    outside = c("nearest", "error"),
    n_threads = NULL
) {
    outside = match.arg(outside)
    if (is.null(segmentation_fit$basis_weights)) {
        stop("segmentation_fit must contain basis_weights.", call. = FALSE)
    }
    if (is.null(segmentation_fit$basis_lattice) || is.null(segmentation_fit$basis_points)) {
        stop("segmentation_fit must contain basis_lattice and basis_points.", call. = FALSE)
    }

    parent_lattice = as.matrix(segmentation_fit$basis_lattice)
    parent_points = as.matrix(segmentation_fit$basis_points)
    basis_weights = as.matrix(segmentation_fit$basis_weights)
    storage.mode(parent_lattice) = "integer"
    storage.mode(parent_points) = "double"
    storage.mode(basis_weights) = "double"

    if (nrow(basis_weights) != nrow(parent_lattice)) {
        stop("segmentation_fit$basis_weights must have one row per parent basis point.", call. = FALSE)
    }
    if (is.null(basis)) {
        basis = if (!is.null(segmentation_fit$parameters$basis)) {
            segmentation_fit$parameters$basis
        } else if (ncol(parent_points) == 2L) {
            "2d"
        } else if (ncol(parent_points) == 3L) {
            "3d"
        } else {
            NA_character_
        }
    }
    basis = normalize_spatial_basis(basis)
    d = if (basis == "2d") 2L else 3L
    if (ncol(parent_lattice) != d || ncol(parent_points) != d) {
        stop("Parent basis geometry does not match selected basis dimensionality.", call. = FALSE)
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
            infer_spatial_basis_origin(parent_lattice, parent_points, basis, s)
        }
    }
    if (length(origin) != d || any(!is.finite(origin))) {
        stop("origin must be a finite numeric vector with length matching the selected basis.", call. = FALSE)
    }
    if (is.null(n_threads)) {
        n_threads = 0L
    } else if (length(n_threads) != 1L || !is.finite(n_threads) || n_threads < 1) {
        stop("n_threads must be NULL or a positive integer.", call. = FALSE)
    } else {
        n_threads = as.integer(n_threads)
    }
    basis_n_threads = if (n_threads == 0L) NULL else n_threads

    coords = as.matrix(coords)
    storage.mode(coords) = "double"
    if (ncol(coords) != d || any(!is.finite(coords))) {
        stop("coords must be a finite matrix with one column per spatial dimension.", call. = FALSE)
    }

    nearest_parent = RANN::nn2(parent_points, query = coords, k = 1L)
    exact_parent = nearest_parent$nn.dists[, 1L] <= sqrt(.Machine$double.eps)
    logits = matrix(NA_real_, nrow = nrow(coords), ncol = ncol(basis_weights))
    if (any(exact_parent)) {
        logits[exact_parent, ] = basis_weights[nearest_parent$nn.idx[exact_parent, 1L], , drop = FALSE]
    }

    bary = if (basis == "2d") {
        tri_barycentric(coords, s = s, origin = origin, n_threads = basis_n_threads)
    } else {
        bcc_barycentric(coords, s = s, origin = origin, n_threads = basis_n_threads)
    }
    design = basis_design_from_barycentric(bary)
    parent_key = do.call(paste, c(as.data.frame(parent_lattice), sep = ":"))
    design_key = do.call(paste, c(as.data.frame(design$basis_lattice), sep = ":"))
    parent_lookup = seq_along(parent_key)
    names(parent_lookup) = parent_key
    matched_vertex = unname(parent_lookup[design_key])
    matched_basis_id = matrix(matched_vertex[design$basis_id], nrow = nrow(design$basis_id))
    row_has_parent_support = rowSums(is.na(matched_basis_id)) == 0L
    interpolate_rows = which(!exact_parent & row_has_parent_support)
    if (length(interpolate_rows) > 0L) {
        logits[interpolate_rows, ] = 0
        for (a in seq_len(ncol(design$basis_id))) {
            logits[interpolate_rows, ] =
                logits[interpolate_rows, ] +
                design$basis_weight[interpolate_rows, a] *
                    basis_weights[matched_basis_id[interpolate_rows, a], , drop = FALSE]
        }
    }
    missing_rows = which(!is.finite(rowSums(logits)))
    if (length(missing_rows) > 0L) {
        if (outside == "error") {
            stop(
                "Could not evaluate the spatial prior at ",
                length(missing_rows),
                " point(s) because their active parent-basis vertices are outside ",
                "segmentation_fit$basis_lattice.",
                call. = FALSE
            )
        }
        logits[missing_rows, ] = basis_weights[nearest_parent$nn.idx[missing_rows, 1L], , drop = FALSE]
    }
    colnames(logits) = colnames(basis_weights)
    prior = softmax_rows(logits)
    colnames(prior) = colnames(basis_weights)
    prior
}

summarize_basis_basins <- function(basin, parent, type_density, cell_types) {
    out = vector("list", length(cell_types))
    for (k in seq_along(cell_types)) {
        basins = sort(unique(basin[, k]))
        idx = match(basin[, k], basins)
        n_basis_points = tabulate(idx, nbins = length(basins))
        mode_density = as.numeric(tapply(type_density[, k], basin[, k], max)[as.character(basins)])
        mode_basis_index = vapply(basins, function(b) {
            which_b = which(basin[, k] == b)
            which_b[which.max(type_density[which_b, k])]
        }, integer(1L))
        out[[k]] = data.frame(
            cell_type = cell_types[k],
            basin = basins,
            n_basis_points = n_basis_points,
            mode_basis_index = mode_basis_index,
            mode_density = mode_density
        )
    }
    do.call(rbind, out)
}

assign_transcripts_to_basis_basins <- function(transcript_basis_id, type_density, basin, cell_types) {
    transcript_basis_id = as.matrix(transcript_basis_id)
    storage.mode(transcript_basis_id) = "integer"
    n = nrow(transcript_basis_id)
    K = ncol(type_density)
    transcript_mode_basis = matrix(NA_integer_, nrow = n, ncol = K, dimnames = list(NULL, cell_types))
    transcript_basin = matrix(NA_integer_, nrow = n, ncol = K, dimnames = list(NULL, cell_types))

    for (k in seq_len(K)) {
        local_density = matrix(type_density[transcript_basis_id, k], nrow = n)
        best_active = max.col(local_density, ties.method = "first")
        mode_basis = transcript_basis_id[cbind(seq_len(n), best_active)]
        transcript_mode_basis[, k] = mode_basis
        transcript_basin[, k] = basin[cbind(mode_basis, k)]
    }

    list(
        transcript_mode_basis = transcript_mode_basis,
        transcript_basin = transcript_basin
    )
}

#' Initial basis-level watershed basins by cell type
#'
#' Runs density ascent on a spatial basis graph separately for each cell type.
#' The per-type watershed score is the product of a fitted total density field
#' and the spatial prior from [spatial_basis_segmentation()]. This function
#' returns only the initial basins; basin merging is intentionally left to a
#' later step.
#'
#' @param segmentation_fit Result from [spatial_basis_segmentation()].
#' @param density_fit Result from [spatial_basis_total_density_field()].
#' @param distance_weight Non-negative penalty for long uphill ascent edges on
#'   the basis graph.
#' @param prior_outside Character; how to evaluate the spatial prior for
#'   density-basis vertices whose parent-basis interpolation would require
#'   parent vertices outside `segmentation_fit$basis_lattice`. `"nearest"` uses
#'   the nearest parent basis vertex logits; `"error"` stops.
#' @param n_threads Integer number of OpenMP threads for prior interpolation.
#'   If `NULL`, uses runtime default.
#' @param show_progress Logical; print progress messages.
#'
#' @return A list with density-basis points, basis edges, interpolated spatial
#'   prior on the density basis, per-cell-type density scores, ascent parents,
#'   initial basin labels, transcript-level basin labels, basin summaries, and
#'   effective parameters.
#' @export
partition_basis_watershed_initial <- function(
    segmentation_fit,
    density_fit,
    distance_weight = 0,
    prior_outside = c("nearest", "error"),
    n_threads = NULL,
    show_progress = TRUE
) {
    prior_outside = match.arg(prior_outside)
    if (is.null(density_fit$total_density_basis)) {
        stop("density_fit must contain total_density_basis.", call. = FALSE)
    }
    if (is.null(density_fit$basis_points) || is.null(density_fit$basis_edges)) {
        stop("density_fit must contain basis_points and basis_edges.", call. = FALSE)
    }
    check_finite_scalar(distance_weight, "distance_weight", lower = 0)

    basis_points = as.matrix(density_fit$basis_points)
    storage.mode(basis_points) = "double"
    total_density = as.numeric(density_fit$total_density_basis)
    if (length(total_density) != nrow(basis_points) || any(!is.finite(total_density)) || any(total_density <= 0)) {
        stop("density_fit$total_density_basis must contain one positive finite value per density basis point.", call. = FALSE)
    }
    basis_edges = as.data.frame(density_fit$basis_edges)
    if (!all(c("from", "to", "distance") %in% colnames(basis_edges))) {
        stop("density_fit$basis_edges must contain from, to, and distance columns.", call. = FALSE)
    }
    if (nrow(basis_edges) == 0L) {
        stop("density_fit$basis_edges must contain at least one edge.", call. = FALSE)
    }

    if (show_progress) {
        message("Evaluating spatial prior on density basis...")
    }
    spatial_prior = evaluate_spatial_prior_on_points(
        segmentation_fit = segmentation_fit,
        coords = basis_points,
        outside = prior_outside,
        n_threads = n_threads
    )
    cell_types = colnames(spatial_prior)
    if (is.null(cell_types)) {
        cell_types = paste0("type", seq_len(ncol(spatial_prior)))
        colnames(spatial_prior) = cell_types
    }

    type_density = sweep(spatial_prior, 1L, total_density, "*")
    colnames(type_density) = cell_types

    if (show_progress) {
        message("Running basis watershed by cell type...")
    }
    n_basis = nrow(basis_points)
    K = ncol(type_density)
    parent = matrix(NA_integer_, nrow = n_basis, ncol = K, dimnames = list(NULL, cell_types))
    basin = matrix(NA_integer_, nrow = n_basis, ncol = K, dimnames = list(NULL, cell_types))
    zero_js = rep(0, nrow(basis_edges))
    for (k in seq_len(K)) {
        ascent = density_ascent_partition_cpp(
            as.integer(basis_edges$from),
            as.integer(basis_edges$to),
            as.numeric(basis_edges$distance),
            type_density[, k],
            zero_js,
            distance_weight,
            0
        )
        parent[, k] = ascent$parent
        basin[, k] = ascent$basin
    }
    transcript_assignment = NULL
    if (!is.null(density_fit$transcript_basis_id)) {
        transcript_assignment = assign_transcripts_to_basis_basins(
            transcript_basis_id = density_fit$transcript_basis_id,
            type_density = type_density,
            basin = basin,
            cell_types = cell_types
        )
    }

    list(
        basis_points = basis_points,
        basis_lattice = density_fit$basis_lattice,
        basis_edges = basis_edges,
        total_density_basis = total_density,
        spatial_prior_basis = spatial_prior,
        type_density_basis = type_density,
        parent = parent,
        basin = basin,
        transcript_mode_basis = transcript_assignment$transcript_mode_basis,
        transcript_basin = transcript_assignment$transcript_basin,
        basin_summary = summarize_basis_basins(basin, parent, type_density, cell_types),
        parameters = list(
            distance_weight = distance_weight,
            prior_outside = prior_outside,
            density_basis_subdivision = density_fit$parameters$basis_subdivision,
            density_quadrature_subdivision = density_fit$parameters$quadrature_subdivision
        )
    )
}
