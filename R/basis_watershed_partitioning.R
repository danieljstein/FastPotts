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

compute_basis_type_support <- function(transcript_basis_id, transcript_basis_weight, posterior, n_basis) {
    transcript_basis_id = as.matrix(transcript_basis_id)
    transcript_basis_weight = as.matrix(transcript_basis_weight)
    posterior = as.matrix(posterior)
    storage.mode(transcript_basis_id) = "integer"
    storage.mode(transcript_basis_weight) = "double"
    storage.mode(posterior) = "double"
    if (!all(dim(transcript_basis_id) == dim(transcript_basis_weight))) {
        stop("transcript_basis_id and transcript_basis_weight must have the same dimensions.", call. = FALSE)
    }
    if (nrow(posterior) != nrow(transcript_basis_id)) {
        stop("posterior must have one row per transcript.", call. = FALSE)
    }

    support = matrix(0, nrow = n_basis, ncol = ncol(posterior))
    for (a in seq_len(ncol(transcript_basis_id))) {
        ok = !is.na(transcript_basis_id[, a]) & transcript_basis_id[, a] >= 1L & transcript_basis_id[, a] <= n_basis
        if (!any(ok)) {
            next
        }
        partial_rows = rowsum(
            posterior[ok, , drop = FALSE] * transcript_basis_weight[ok, a],
            group = transcript_basis_id[ok, a],
            reorder = FALSE
        )
        support[as.integer(rownames(partial_rows)), ] = support[as.integer(rownames(partial_rows)), ] + partial_rows
    }
    colnames(support) = colnames(posterior)
    support
}

compute_basis_type_support_from_segmentation_fit <- function(segmentation_fit, density_fit, cell_types, n_threads = NULL) {
    if (is.null(segmentation_fit$transcripts_df)) {
        stop(
            "posterior is NULL and segmentation_fit does not contain transcripts_df; ",
            "provide posterior explicitly or keep transcripts_df in the segmentation fit.",
            call. = FALSE
        )
    }
    if (is.null(segmentation_fit$basis_weights) || is.null(segmentation_fit$basis_lattice)) {
        stop("segmentation_fit must contain basis_weights and basis_lattice.", call. = FALSE)
    }
    if (is.null(segmentation_fit$cell_signatures)) {
        stop("segmentation_fit must contain cell_signatures.", call. = FALSE)
    }

    transcripts_df = segmentation_fit$transcripts_df
    n = nrow(density_fit$transcript_basis_id)
    if (nrow(transcripts_df) != n) {
        stop(
            "posterior is NULL and segmentation_fit$transcripts_df does not have one row per density-fit transcript; ",
            "provide a posterior matrix explicitly.",
            call. = FALSE
        )
    }

    basis = normalize_spatial_basis(segmentation_fit$parameters$basis)
    x = segmentation_fit$parameters$x %||% "x_location"
    y = segmentation_fit$parameters$y %||% "y_location"
    z = segmentation_fit$parameters$z %||% "z_location"
    gene = segmentation_fit$parameters$gene %||% "feature_name"
    coord_cols = if (basis == "2d") c(x, y) else c(x, y, z)
    missing_cols = setdiff(c(coord_cols, gene), colnames(transcripts_df))
    if (length(missing_cols) > 0L) {
        stop("Missing required column(s) in segmentation_fit$transcripts_df: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    gene_index = match(transcripts_df[[gene]], rownames(segmentation_fit$cell_signatures))
    if (anyNA(gene_index)) {
        stop("Some transcript genes are absent from segmentation_fit$cell_signatures.", call. = FALSE)
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
    bary = if (basis == "2d") {
        tri_barycentric(
            coords,
            s = segmentation_fit$parameters$s,
            origin = segmentation_fit$parameters$origin,
            n_threads = basis_n_threads
        )
    } else {
        bcc_barycentric(
            coords,
            s = segmentation_fit$parameters$s,
            origin = segmentation_fit$parameters$origin,
            n_threads = basis_n_threads
        )
    }
    design = remap_spatial_basis_design(basis_design_from_barycentric(bary), segmentation_fit$basis_lattice)

    support = spatial_basis_density_support_cpp(
        par = as.vector(as.matrix(segmentation_fit$basis_weights)),
        segmentation_basis_id = design$basis_id - 1L,
        segmentation_basis_weight = design$basis_weight,
        gene_index = as.integer(gene_index - 1L),
        log_signature = log(segmentation_fit$cell_signatures),
        density_basis_id = as.matrix(density_fit$transcript_basis_id) - 1L,
        density_basis_weight = as.matrix(density_fit$transcript_basis_weight),
        n_segmentation_basis = nrow(segmentation_fit$basis_weights),
        n_density_basis = nrow(density_fit$basis_points),
        n_threads = n_threads,
        n_cell_types = ncol(segmentation_fit$cell_signatures)
    )
    colnames(support) = colnames(segmentation_fit$cell_signatures)
    support[, cell_types, drop = FALSE]
}

summarize_active_basis_basins <- function(basin, mode_basis, active_start, type_density, cell_types) {
    out = vector("list", length(cell_types))
    for (k in seq_along(cell_types)) {
        cell_type = cell_types[k]
        basins = sort(unique(stats::na.omit(basin[[cell_type]])))
        if (length(basins) == 0L) {
            out[[k]] = NULL
            next
        }
        n_basis_points = tabulate(match(basin[[cell_type]][active_start[[cell_type]]], basins), nbins = length(basins))
        mode = mode_basis[[cell_type]][basins]
        out[[k]] = data.frame(
            cell_type = cell_type,
            basin = basins,
            cell = paste(cell_type, mode, sep = "-"),
            n_active_basis_points = n_basis_points,
            mode_basis_index = mode,
            mode_density = type_density[[cell_type]][mode]
        )
    }
    do.call(rbind, out)
}

make_basis_cell_names <- function(cell_type, basin, mode_basis) {
    out = rep(NA_character_, length(basin))
    ok = !is.na(basin)
    out[ok] = paste(cell_type, mode_basis[basin[ok]], sep = "-")
    out
}

assign_transcripts_to_basis_basins <- function(transcript_basis_id, type_density, basin, active_start, cell_types) {
    transcript_basis_id = as.matrix(transcript_basis_id)
    storage.mode(transcript_basis_id) = "integer"
    n = nrow(transcript_basis_id)
    K = length(cell_types)
    transcript_mode_basis = matrix(NA_integer_, nrow = n, ncol = K, dimnames = list(NULL, cell_types))
    transcript_basin = matrix(NA_integer_, nrow = n, ncol = K, dimnames = list(NULL, cell_types))

    for (k in seq_len(K)) {
        local_basis_id = transcript_basis_id
        local_basis_id[is.na(local_basis_id)] = NA_integer_
        local_density = matrix(type_density[[k]][local_basis_id], nrow = n)
        local_active = matrix(active_start[[k]][local_basis_id], nrow = n)
        local_active[is.na(local_active)] = FALSE
        local_density[!local_active] = -Inf
        has_active = rowSums(local_active) > 0L
        best_active = max.col(local_density[has_active, , drop = FALSE], ties.method = "first")
        mode_basis = transcript_basis_id[has_active, , drop = FALSE][cbind(seq_len(sum(has_active)), best_active)]
        transcript_mode_basis[has_active, k] = mode_basis
        transcript_basin[has_active, k] = basin[[k]][mode_basis]
    }

    list(
        transcript_mode_basis = transcript_mode_basis,
        transcript_basin = transcript_basin
    )
}

basin_cell_lookup <- function(basis_partition) {
    summary = basis_partition$basin_summary
    out = vector("list", length(basis_partition$active_cell_types))
    names(out) = basis_partition$active_cell_types
    for (cell_type in basis_partition$active_cell_types) {
        rows = summary$cell_type == cell_type
        values = as.character(summary$cell[rows])
        names(values) = as.character(summary$basin[rows])
        out[[cell_type]] = values
    }
    out
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
#' @param posterior Optional transcript-by-cell-type posterior matrix. If
#'   `NULL`, uses `segmentation_fit$marginals` when present; otherwise streams
#'   posterior support from compact `segmentation_fit` fields without
#'   materializing the full posterior matrix.
#' @param min_posterior_support Minimum posterior-weighted transcript support
#'   for a density-basis vertex to be an active watershed start for a cell type.
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
    posterior = NULL,
    min_posterior_support = 0,
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
    check_finite_scalar(min_posterior_support, "min_posterior_support", lower = 0)
    if (is.null(density_fit$transcript_basis_id) || is.null(density_fit$transcript_basis_weight)) {
        stop("density_fit must contain transcript_basis_id and transcript_basis_weight.", call. = FALSE)
    }

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

    if (is.null(posterior) && !is.null(segmentation_fit$marginals)) {
        posterior = segmentation_fit$marginals
    }
    if (is.null(posterior)) {
        if (show_progress) {
            message("Computing posterior support on density basis...")
        }
        support = compute_basis_type_support_from_segmentation_fit(
            segmentation_fit = segmentation_fit,
            density_fit = density_fit,
            cell_types = cell_types,
            n_threads = n_threads
        )
    } else {
        posterior = normalize_posterior_matrix(posterior, nrow(density_fit$transcript_basis_id))
        if (is.null(colnames(posterior))) {
            colnames(posterior) = cell_types
        }
        if (!identical(colnames(posterior), cell_types)) {
            posterior = posterior[, cell_types, drop = FALSE]
        }
        support = compute_basis_type_support(
            transcript_basis_id = density_fit$transcript_basis_id,
            transcript_basis_weight = density_fit$transcript_basis_weight,
            posterior = posterior,
            n_basis = nrow(basis_points)
        )
    }
    active_start_matrix = support > min_posterior_support
    active_cell_types = cell_types[colSums(active_start_matrix) > 0L]
    if (length(active_cell_types) == 0L) {
        stop("No cell types had active basis starts at the requested min_posterior_support.", call. = FALSE)
    }

    if (show_progress) {
        message("Running basis watershed by cell type...")
    }
    n_basis = nrow(basis_points)
    parent = vector("list", length(active_cell_types))
    root = vector("list", length(active_cell_types))
    basin = vector("list", length(active_cell_types))
    mode_basis = vector("list", length(active_cell_types))
    type_density = vector("list", length(active_cell_types))
    active_start = vector("list", length(active_cell_types))
    names(parent) = names(root) = names(basin) = names(mode_basis) =
        names(type_density) = names(active_start) = active_cell_types
    for (cell_type in active_cell_types) {
        type_density[[cell_type]] = total_density * spatial_prior[, cell_type]
        active_start[[cell_type]] = active_start_matrix[, cell_type]
        ascent = density_ascent_active_partition_cpp(
            as.integer(basis_edges$from),
            as.integer(basis_edges$to),
            as.numeric(basis_edges$distance),
            type_density[[cell_type]],
            active_start[[cell_type]],
            distance_weight
        )
        parent[[cell_type]] = ascent$parent
        root[[cell_type]] = ascent$root
        basin[[cell_type]] = ascent$basin
        mode_basis[[cell_type]] = ascent$mode_basis
    }
    transcript_assignment = assign_transcripts_to_basis_basins(
        transcript_basis_id = density_fit$transcript_basis_id,
        type_density = type_density,
        basin = basin,
        active_start = active_start,
        cell_types = active_cell_types
    )

    list(
        basis_points = basis_points,
        basis_lattice = density_fit$basis_lattice,
        basis_edges = basis_edges,
        total_density_basis = total_density,
        spatial_prior_basis = spatial_prior,
        posterior_support_basis = support,
        active_start = active_start,
        active_cell_types = active_cell_types,
        type_density_basis = type_density,
        parent = parent,
        root = root,
        basin = basin,
        mode_basis = mode_basis,
        transcript_mode_basis = transcript_assignment$transcript_mode_basis,
        transcript_basin = transcript_assignment$transcript_basin,
        basin_summary = summarize_active_basis_basins(basin, mode_basis, active_start, type_density, active_cell_types),
        parameters = list(
            min_posterior_support = min_posterior_support,
            distance_weight = distance_weight,
            prior_outside = prior_outside,
            density_basis_subdivision = density_fit$parameters$basis_subdivision,
            density_quadrature_subdivision = density_fit$parameters$quadrature_subdivision
        )
    )
}

#' Initial basis-level watershed basins from total density only
#'
#' Runs density ascent on a spatial basis graph using only
#' `density_fit$total_density_basis`, without a segmentation fit, spatial prior,
#' or transcript cell-type posteriors. This is useful when the desired
#' partition should follow total transcript density rather than cell-type-
#' specific density.
#'
#' @param density_fit Result from [spatial_basis_log_density_field()] or
#'   [spatial_basis_total_density_field()] containing `total_density_basis`,
#'   `basis_points`, `basis_edges`, and transcript basis interpolation.
#' @param distance_weight Non-negative penalty for long uphill ascent edges on
#'   the basis graph.
#' @param cell_type Character label used for the single total-density watershed
#'   column in transcript assignments and count matrices.
#' @param show_progress Logical; print progress messages.
#'
#' @return A basis-watershed partition object compatible with
#'   [basis_watershed_gene_counts()]. It contains one active cell type named by
#'   `cell_type`.
#' @export
partition_basis_watershed_total <- function(
    density_fit,
    distance_weight = 0,
    cell_type = "total",
    show_progress = TRUE
) {
    if (is.null(density_fit$total_density_basis)) {
        stop("density_fit must contain total_density_basis.", call. = FALSE)
    }
    if (is.null(density_fit$basis_points) || is.null(density_fit$basis_edges)) {
        stop("density_fit must contain basis_points and basis_edges.", call. = FALSE)
    }
    if (is.null(density_fit$transcript_basis_id)) {
        stop("density_fit must contain transcript_basis_id.", call. = FALSE)
    }
    check_finite_scalar(distance_weight, "distance_weight", lower = 0)
    if (length(cell_type) != 1L || is.na(cell_type) || !nzchar(cell_type)) {
        stop("cell_type must be a non-empty character scalar.", call. = FALSE)
    }

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
        message("Running basis watershed on total density...")
    }
    n_basis = nrow(basis_points)
    active = rep(TRUE, n_basis)
    ascent = density_ascent_active_partition_cpp(
        as.integer(basis_edges$from),
        as.integer(basis_edges$to),
        as.numeric(basis_edges$distance),
        total_density,
        active,
        distance_weight
    )

    cell_types = cell_type
    parent = root = basin = mode_basis = type_density = active_start = list()
    parent[[cell_type]] = ascent$parent
    root[[cell_type]] = ascent$root
    basin[[cell_type]] = ascent$basin
    mode_basis[[cell_type]] = ascent$mode_basis
    type_density[[cell_type]] = total_density
    active_start[[cell_type]] = active

    transcript_assignment = assign_transcripts_to_basis_basins(
        transcript_basis_id = density_fit$transcript_basis_id,
        type_density = type_density,
        basin = basin,
        active_start = active_start,
        cell_types = cell_types
    )

    list(
        basis_points = basis_points,
        basis_lattice = density_fit$basis_lattice,
        basis_edges = basis_edges,
        total_density_basis = total_density,
        spatial_prior_basis = NULL,
        posterior_support_basis = NULL,
        active_start = active_start,
        active_cell_types = cell_types,
        type_density_basis = type_density,
        parent = parent,
        root = root,
        basin = basin,
        mode_basis = mode_basis,
        transcript_mode_basis = transcript_assignment$transcript_mode_basis,
        transcript_basin = transcript_assignment$transcript_basin,
        basin_summary = summarize_active_basis_basins(basin, mode_basis, active_start, type_density, cell_types),
        parameters = list(
            mode = "total_density",
            distance_weight = distance_weight,
            density_basis_subdivision = density_fit$parameters$basis_subdivision,
            density_quadrature_subdivision = density_fit$parameters$quadrature_subdivision
        )
    )
}

#' Collect transcript genes into basis-watershed basin counts
#'
#' Builds a sparse gene-by-basin count matrix from the transcript assignments
#' returned by [partition_basis_watershed_initial()].
#'
#' @param basis_partition Result from [partition_basis_watershed_initial()].
#' @param transcripts_df Transcript-level data frame corresponding to the
#'   density and segmentation fits.
#' @param posterior Optional transcript-by-cell-type posterior matrix. Required
#'   for `mode = "weighted"` and used for max-posterior assignment in
#'   `mode = "max_posterior"`. If `NULL`, `mode = "max_posterior"` uses the
#'   active cell type with the largest basis-watershed type density at each
#'   transcript's selected mode basis point.
#' @param gene Character gene column name.
#' @param mode Character; `"max_posterior"` assigns each transcript once to its
#'   maximum-posterior cell type, while `"weighted"` contributes posterior
#'   weight to each active cell type with a valid transcript basin.
#'
#' @return Sparse `dgCMatrix` with genes in rows and `celltype-modebasis`
#'   initial watershed cells in columns.
#' @export
basis_watershed_gene_counts <- function(
    basis_partition,
    transcripts_df,
    posterior = NULL,
    gene = "feature_name",
    mode = c("max_posterior", "weighted")
) {
    mode = match.arg(mode)
    if (!gene %in% colnames(transcripts_df)) {
        stop("Column '", gene, "' not found in transcripts_df.", call. = FALSE)
    }
    if (is.null(basis_partition$transcript_basin) || is.null(basis_partition$basin_summary)) {
        stop("basis_partition must contain transcript_basin and basin_summary.", call. = FALSE)
    }
    cell_types = colnames(basis_partition$transcript_basin)
    n = nrow(basis_partition$transcript_basin)
    if (nrow(transcripts_df) != n) {
        stop("transcripts_df must have one row per transcript assignment.", call. = FALSE)
    }
    total_only = mode == "max_posterior" &&
        length(cell_types) == 1L &&
        isTRUE(basis_partition$parameters$mode == "total_density")
    if (!is.null(posterior) && mode == "weighted") {
        posterior = normalize_posterior_matrix(posterior, n)
        if (is.null(colnames(posterior))) {
            stop("posterior must have column names matching cell types.", call. = FALSE)
        }
    } else if (!is.null(posterior) && !total_only) {
        posterior = as.matrix(posterior)
        storage.mode(posterior) = "double"
        if (nrow(posterior) != n) {
            stop("posterior must have one row per transcript.", call. = FALSE)
        }
        if (ncol(posterior) < 1L) {
            stop("posterior must have at least one column.", call. = FALSE)
        }
        if (any(!is.finite(posterior)) || any(posterior < 0)) {
            stop("posterior must contain finite non-negative values.", call. = FALSE)
        }
        if (is.null(colnames(posterior))) {
            stop("posterior must have column names matching cell types.", call. = FALSE)
        }
    } else if (mode == "weighted") {
        stop("posterior is required when mode = 'weighted'.", call. = FALSE)
    }

    gene_values = as.character(transcripts_df[[gene]])
    gene_levels = sort(unique(gene_values))
    gene_id = match(gene_values, gene_levels)
    cell_levels = as.character(basis_partition$basin_summary$cell)
    cell_id_lookup = seq_along(cell_levels)
    names(cell_id_lookup) = cell_levels
    basin_lookup = basin_cell_lookup(basis_partition)

    ii = list()
    jj = list()
    xx = list()
    part_i = 0L

    if (mode == "max_posterior") {
        if (total_only) {
            max_type = rep(cell_types, n)
        } else if (!is.null(posterior)) {
            max_type = colnames(posterior)[max.col(posterior, ties.method = "first")]
        } else {
            score = matrix(-Inf, nrow = n, ncol = length(cell_types), dimnames = list(NULL, cell_types))
            for (cell_type in cell_types) {
                mode_basis = basis_partition$transcript_mode_basis[, cell_type]
                ok = !is.na(mode_basis)
                score[ok, cell_type] = basis_partition$type_density_basis[[cell_type]][mode_basis[ok]]
            }
            max_type = colnames(score)[max.col(score, ties.method = "first")]
        }
        for (cell_type in cell_types) {
            rows = which(max_type == cell_type & !is.na(basis_partition$transcript_basin[, cell_type]))
            if (length(rows) == 0L) next
            cell_name = basin_lookup[[cell_type]][as.character(basis_partition$transcript_basin[rows, cell_type])]
            keep = !is.na(cell_name)
            rows = rows[keep]
            cell_name = cell_name[keep]
            part_i = part_i + 1L
            ii[[part_i]] = gene_id[rows]
            jj[[part_i]] = unname(cell_id_lookup[cell_name])
            xx[[part_i]] = rep(1, length(rows))
        }
    } else {
        for (cell_type in cell_types) {
            if (!cell_type %in% colnames(posterior)) next
            rows = which(!is.na(basis_partition$transcript_basin[, cell_type]) & posterior[, cell_type] > 0)
            if (length(rows) == 0L) next
            cell_name = basin_lookup[[cell_type]][as.character(basis_partition$transcript_basin[rows, cell_type])]
            keep = !is.na(cell_name)
            rows = rows[keep]
            cell_name = cell_name[keep]
            part_i = part_i + 1L
            ii[[part_i]] = gene_id[rows]
            jj[[part_i]] = unname(cell_id_lookup[cell_name])
            xx[[part_i]] = posterior[rows, cell_type]
        }
    }
    ii = if (length(ii) > 0L) unlist(ii, use.names = FALSE) else integer()
    jj = if (length(jj) > 0L) unlist(jj, use.names = FALSE) else integer()
    xx = if (length(xx) > 0L) unlist(xx, use.names = FALSE) else numeric()

    Matrix::sparseMatrix(
        i = ii,
        j = jj,
        x = xx,
        dims = c(length(gene_levels), length(cell_levels)),
        dimnames = list(gene_levels, cell_levels)
    )
}
