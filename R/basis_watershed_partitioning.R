softmax_rows <- function(x) {
    x = as.matrix(x)
    row_max = do.call(pmax, as.data.frame(x))
    z = exp(sweep(x, 1L, row_max, "-"))
    sweep(z, 1L, rowSums(z), "/")
}

normalize_cell_type_groups <- function(cell_type_groups, fine_cell_types) {
    if (is.null(cell_type_groups)) {
        out = as.list(fine_cell_types)
        names(out) = fine_cell_types
        return(out)
    }
    if (!is.list(cell_type_groups) || is.null(names(cell_type_groups)) || any(!nzchar(names(cell_type_groups)))) {
        stop("cell_type_groups must be a named list of fine cell-type names.", call. = FALSE)
    }
    if (anyDuplicated(names(cell_type_groups))) {
        stop("cell_type_groups names must be unique.", call. = FALSE)
    }
    out = lapply(cell_type_groups, as.character)
    empty = vapply(out, length, integer(1L)) == 0L
    if (any(empty)) {
        stop("Every cell_type_groups entry must contain at least one fine cell type.", call. = FALSE)
    }
    members = unlist(out, use.names = FALSE)
    missing = setdiff(members, fine_cell_types)
    if (length(missing) > 0L) {
        stop("cell_type_groups contains unknown fine cell type(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    duplicate_members = unique(members[duplicated(members)])
    if (length(duplicate_members) > 0L) {
        stop("Fine cell type(s) appear in more than one cell_type_groups entry: ", paste(duplicate_members, collapse = ", "), call. = FALSE)
    }
    unmapped = setdiff(fine_cell_types, members)
    if (length(unmapped) > 0L) {
        singleton = as.list(unmapped)
        names(singleton) = unmapped
        out = c(out, singleton)
    }
    if (anyDuplicated(names(out))) {
        stop("Coarse cell type names must be unique and must not duplicate unmapped fine cell type names.", call. = FALSE)
    }
    out
}

aggregate_matrix_by_cell_type_groups <- function(x, cell_type_groups) {
    x = as.matrix(x)
    if (is.null(colnames(x))) {
        stop("Matrix to aggregate must have column names.", call. = FALSE)
    }
    out = matrix(0, nrow = nrow(x), ncol = length(cell_type_groups))
    colnames(out) = names(cell_type_groups)
    for (k in seq_along(cell_type_groups)) {
        cols = cell_type_groups[[k]]
        if (length(cols) == 1L) {
            out[, k] = x[, cols]
        } else {
            out[, k] = rowSums(x[, cols, drop = FALSE])
        }
    }
    out
}

prepare_grouped_posterior <- function(posterior, n, target_cell_types, cell_type_groups = NULL) {
    posterior = normalize_posterior_matrix(posterior, n)
    if (is.null(colnames(posterior))) {
        if (is.null(cell_type_groups)) {
            colnames(posterior) = target_cell_types
        } else {
            stop("posterior must have column names when cell type grouping is used.", call. = FALSE)
        }
    }
    if (all(target_cell_types %in% colnames(posterior))) {
        return(posterior[, target_cell_types, drop = FALSE])
    }
    if (!is.null(cell_type_groups)) {
        groups = normalize_cell_type_groups(cell_type_groups, colnames(posterior))
        posterior = aggregate_matrix_by_cell_type_groups(posterior, groups)
    }
    missing = setdiff(target_cell_types, colnames(posterior))
    if (length(missing) > 0L) {
        stop("posterior is missing required cell type(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    posterior[, target_cell_types, drop = FALSE]
}

partition_cell_type_groups <- function(basis_partition) {
    basis_partition$parameters$cell_type_groups %||% NULL
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

#' Basis-level spatial overlap between segmentation signatures
#'
#' Computes pairwise overlap between cell-signature spatial prior fields from a
#' fitted [spatial_basis_segmentation()] object. When `density_fit` is supplied,
#' overlap is evaluated on the density basis and weighted by
#' `density_fit$total_density_basis`; otherwise it is evaluated on the
#' segmentation basis with uniform basis-point weights.
#'
#' This diagnostic is intended as a cheap pre-partition screen for fine
#' signatures that occupy the same spatial support and may be better grouped
#' into a coarser cell type before basis watershed partitioning.
#'
#' @param segmentation_fit Result from [spatial_basis_segmentation()].
#' @param density_fit Optional result from [spatial_basis_total_density_field()]
#'   or a compatible object containing `basis_points` and `total_density_basis`.
#' @param cell_type_groups Optional named list mapping coarse cell type names to
#'   fine cell-signature names. Fine signatures not listed in any group are
#'   retained as singleton groups before overlap is computed.
#' @param cell_types Optional character vector of cell types to include after
#'   grouping. Defaults to all available cell types.
#' @param prior_outside Character; how to evaluate the spatial prior for
#'   density-basis vertices whose parent-basis interpolation would require
#'   parent vertices outside `segmentation_fit$basis_lattice`. `"nearest"` uses
#'   the nearest parent basis vertex logits; `"error"` stops.
#' @param n_threads Integer number of OpenMP threads for prior interpolation.
#'   If `NULL`, uses runtime default.
#'
#' @return A data frame with one row per cell-type pair. `shared_mass` is the
#'   density-weighted sum of the pointwise minimum of the two spatial priors.
#'   `overlap_coef` is `shared_mass / min(mass_a, mass_b)`, so values near 1
#'   mean most of the smaller signature's spatial mass is co-localized with the
#'   other signature. `jaccard` is `shared_mass / (mass_a + mass_b -
#'   shared_mass)`, `cosine` is a weighted cosine similarity between prior
#'   fields, and `frac_a_shared` / `frac_b_shared` are directional shared-mass
#'   fractions.
#' @export
spatial_basis_signature_overlap <- function(
    segmentation_fit,
    density_fit = NULL,
    cell_type_groups = NULL,
    cell_types = NULL,
    prior_outside = c("nearest", "error"),
    n_threads = NULL
) {
    prior_outside = match.arg(prior_outside)
    if (is.null(segmentation_fit$basis_weights) || is.null(segmentation_fit$basis_points)) {
        stop("segmentation_fit must contain basis_weights and basis_points.", call. = FALSE)
    }

    if (is.null(density_fit)) {
        coords = as.matrix(segmentation_fit$basis_points)
        weights = rep(1, nrow(coords))
        basis_source = "segmentation"
    } else {
        if (is.null(density_fit$basis_points) || is.null(density_fit$total_density_basis)) {
            stop("density_fit must contain basis_points and total_density_basis.", call. = FALSE)
        }
        coords = as.matrix(density_fit$basis_points)
        weights = as.numeric(density_fit$total_density_basis)
        basis_source = "density"
    }
    storage.mode(coords) = "double"
    if (any(!is.finite(coords))) {
        stop("basis points must contain only finite coordinates.", call. = FALSE)
    }
    if (
        length(weights) != nrow(coords) ||
        any(!is.finite(weights)) ||
        any(weights < 0) ||
        sum(weights) <= 0
    ) {
        stop("basis weights must contain one non-negative finite value per basis point and have positive total mass.", call. = FALSE)
    }

    spatial_prior = if (is.null(density_fit)) {
        softmax_rows(segmentation_fit$basis_weights)
    } else {
        evaluate_spatial_prior_on_points(
            segmentation_fit = segmentation_fit,
            coords = coords,
            outside = prior_outside,
            n_threads = n_threads
        )
    }
    fine_cell_types = colnames(spatial_prior)
    if (is.null(fine_cell_types)) {
        fine_cell_types = paste0("type", seq_len(ncol(spatial_prior)))
        colnames(spatial_prior) = fine_cell_types
    }
    group_list = normalize_cell_type_groups(cell_type_groups, fine_cell_types)
    if (!is.null(cell_type_groups)) {
        spatial_prior = aggregate_matrix_by_cell_type_groups(spatial_prior, group_list)
    }

    available_cell_types = colnames(spatial_prior)
    if (is.null(cell_types)) {
        cell_types = available_cell_types
    } else {
        cell_types = as.character(cell_types)
        missing = setdiff(cell_types, available_cell_types)
        if (length(missing) > 0L) {
            stop("cell_types contains unknown cell type(s): ", paste(missing, collapse = ", "), call. = FALSE)
        }
    }
    if (length(cell_types) < 2L) {
        stop("At least two cell types are required to compute pairwise overlap.", call. = FALSE)
    }
    spatial_prior = spatial_prior[, cell_types, drop = FALSE]

    overlap = spatial_basis_signature_overlap_cpp(
        spatial_prior = spatial_prior,
        weights = weights,
        n_threads = if (is.null(n_threads)) 0L else as.integer(n_threads)
    )
    mass_a = overlap$mass[overlap$pair_a]
    mass_b = overlap$mass[overlap$pair_b]
    shared_mass = overlap$shared
    union_mass = mass_a + mass_b - shared_mass
    cosine_denom = sqrt(overlap$norm2[overlap$pair_a] * overlap$norm2[overlap$pair_b])
    out = data.frame(
        cell_type_a = cell_types[overlap$pair_a],
        cell_type_b = cell_types[overlap$pair_b],
        mass_a = mass_a,
        mass_b = mass_b,
        shared_mass = shared_mass,
        overlap_coef = shared_mass / pmin(mass_a, mass_b),
        frac_a_shared = shared_mass / mass_a,
        frac_b_shared = shared_mass / mass_b,
        jaccard = shared_mass / union_mass,
        cosine = ifelse(cosine_denom > 0, overlap$cross / cosine_denom, NA_real_),
        stringsAsFactors = FALSE
    )
    out = out[order(out$overlap_coef, out$shared_mass, decreasing = TRUE), , drop = FALSE]
    rownames(out) = NULL
    attr(out, "basis_source") = basis_source
    attr(out, "n_basis_points") = nrow(coords)
    attr(out, "total_weight") = sum(weights)
    attr(out, "cell_type_groups") = if (is.null(cell_type_groups)) NULL else group_list
    out
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

compute_basis_type_support_max_posterior <- function(transcript_basis_id, transcript_basis_weight, posterior, n_basis) {
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
    max_type = max.col(posterior, ties.method = "first")
    for (a in seq_len(ncol(transcript_basis_id))) {
        ok = !is.na(transcript_basis_id[, a]) & transcript_basis_id[, a] >= 1L & transcript_basis_id[, a] <= n_basis
        if (!any(ok)) {
            next
        }
        for (k in seq_len(ncol(posterior))) {
            ok_k = ok & max_type == k
            if (!any(ok_k)) {
                next
            }
            partial_rows = rowsum(
                transcript_basis_weight[ok_k, a],
                group = transcript_basis_id[ok_k, a],
                reorder = FALSE
            )
            support[as.integer(rownames(partial_rows)), k] = support[as.integer(rownames(partial_rows)), k] +
                as.numeric(partial_rows[, 1L])
        }
    }
    colnames(support) = colnames(posterior)
    support
}

compute_basis_type_support_from_segmentation_fit <- function(segmentation_fit, density_fit, cell_types, n_threads = NULL, support_mode = "posterior") {
    support_mode = match.arg(support_mode, c("posterior", "max_posterior"))
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
        n_cell_types = ncol(segmentation_fit$cell_signatures),
        hard_max = identical(support_mode, "max_posterior")
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

validate_basis_watershed_edges <- function(basis_edges, n_basis) {
    if (!all(c("from", "to", "distance") %in% colnames(basis_edges))) {
        stop("density_fit$basis_edges must contain from, to, and distance columns.", call. = FALSE)
    }
    if (nrow(basis_edges) == 0L) {
        stop("density_fit$basis_edges must contain at least one edge.", call. = FALSE)
    }
    edge_from = as.integer(basis_edges$from)
    edge_to = as.integer(basis_edges$to)
    edge_distance = as.numeric(basis_edges$distance)
    bad_edge = is.na(edge_from) | is.na(edge_to) |
        edge_from < 1L | edge_from > n_basis |
        edge_to < 1L | edge_to > n_basis
    if (any(bad_edge)) {
        bad_i = which(bad_edge)[1L]
        stop(
            "density_fit$basis_edges contains an out-of-range edge endpoint at row ",
            bad_i,
            ". Endpoints must be in 1:nrow(density_fit$basis_points) = 1:",
            n_basis,
            ".",
            call. = FALSE
        )
    }
    if (any(!is.finite(edge_distance) | edge_distance <= 0)) {
        stop("density_fit$basis_edges$distance must contain positive finite values.", call. = FALSE)
    }
    list(from = edge_from, to = edge_to, distance = edge_distance)
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
#' @param cell_type_groups Optional named list mapping coarse watershed cell
#'   type names to fine cell-signature names. Fine signatures not listed in any
#'   group are retained as singleton groups. Spatial priors and posteriors are
#'   summed within each group before active starts and max-posterior support are
#'   computed.
#' @param min_posterior_support Minimum posterior-weighted transcript support
#'   for a density-basis vertex to be an active watershed start for a cell type.
#' @param support_mode Character; `"max_posterior"` first assigns each
#'   transcript to its maximum-posterior cell type and contributes support only
#'   to that type, while `"posterior"` uses soft posterior support for active
#'   starts.
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
    cell_type_groups = NULL,
    min_posterior_support = 0,
    support_mode = c("max_posterior", "posterior"),
    distance_weight = 0,
    prior_outside = c("nearest", "error"),
    n_threads = NULL,
    show_progress = TRUE
) {
    support_mode = match.arg(support_mode)
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
    checked_edges = validate_basis_watershed_edges(basis_edges, nrow(basis_points))

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
    fine_cell_types = cell_types
    group_list = normalize_cell_type_groups(cell_type_groups, fine_cell_types)
    if (!is.null(cell_type_groups)) {
        spatial_prior = aggregate_matrix_by_cell_type_groups(spatial_prior, group_list)
        cell_types = colnames(spatial_prior)
    }

    if (is.null(posterior) && !is.null(segmentation_fit$marginals)) {
        posterior = segmentation_fit$marginals
    }
    if (is.null(posterior)) {
        if (!is.null(cell_type_groups) && identical(support_mode, "max_posterior")) {
            stop(
                "cell_type_groups with support_mode = 'max_posterior' requires a posterior matrix ",
                "or segmentation_fit$marginals so fine posteriors can be summed before choosing the max coarse group.",
                call. = FALSE
            )
        }
        if (show_progress) {
            if (identical(support_mode, "max_posterior")) {
                message("Computing max-posterior support on density basis...")
            } else {
                message("Computing posterior support on density basis...")
            }
        }
        support = compute_basis_type_support_from_segmentation_fit(
            segmentation_fit = segmentation_fit,
            density_fit = density_fit,
            cell_types = fine_cell_types,
            n_threads = n_threads,
            support_mode = support_mode
        )
        if (!is.null(cell_type_groups)) {
            support = aggregate_matrix_by_cell_type_groups(support, group_list)
        }
    } else {
        posterior = prepare_grouped_posterior(
            posterior = posterior,
            n = nrow(density_fit$transcript_basis_id),
            target_cell_types = cell_types,
            cell_type_groups = cell_type_groups
        )
        support_fun = if (identical(support_mode, "max_posterior")) {
            compute_basis_type_support_max_posterior
        } else {
            compute_basis_type_support
        }
        support = support_fun(
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
        if (any(!is.finite(type_density[[cell_type]]))) {
            stop("Computed type density contains non-finite values for cell type '", cell_type, "'.", call. = FALSE)
        }
        active_start[[cell_type]] = active_start_matrix[, cell_type]
        ascent = density_ascent_active_partition_cpp(
            checked_edges$from,
            checked_edges$to,
            checked_edges$distance,
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
            support_mode = support_mode,
            distance_weight = distance_weight,
            prior_outside = prior_outside,
            cell_type_groups = if (is.null(cell_type_groups)) NULL else group_list,
            fine_cell_types = fine_cell_types,
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
#' @param density_fit Result from [spatial_basis_log_density_field()],
#'   `spatial_basis_linear_density_field()`, or
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
    checked_edges = validate_basis_watershed_edges(basis_edges, nrow(basis_points))

    if (show_progress) {
        message("Running basis watershed on total density...")
    }
    n_basis = nrow(basis_points)
    active = rep(TRUE, n_basis)
    ascent = density_ascent_active_partition_cpp(
        checked_edges$from,
        checked_edges$to,
        checked_edges$distance,
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
#'   for `mode = "max_posterior"` and `mode = "weighted"`.
#' @param gene Character gene column name.
#' @param mode Character; `"max_posterior"` assigns each transcript once to its
#'   maximum-posterior cell type, `"weighted"` contributes posterior weight to
#'   each active cell type with a valid transcript basin, and `"density"`
#'   explicitly uses the active cell type with the largest basis-watershed type
#'   density at each transcript's selected mode basis point.
#' @param cell_type_groups Optional named list mapping coarse watershed cell
#'   type names to fine posterior column names. Defaults to the mapping stored
#'   in `basis_partition`, if present. Fine posterior columns are summed within
#'   groups before max-posterior or weighted assignment.
#'
#' @return Sparse `dgCMatrix` with genes in rows and `celltype-modebasis`
#'   initial watershed cells in columns.
#' @export
basis_watershed_gene_counts <- function(
    basis_partition,
    transcripts_df,
    posterior = NULL,
    gene = "feature_name",
    mode = c("max_posterior", "weighted", "density"),
    cell_type_groups = partition_cell_type_groups(basis_partition)
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
    total_only = mode == "density" &&
        length(cell_types) == 1L &&
        isTRUE(basis_partition$parameters$mode == "total_density")
    if (is.null(posterior) && mode %in% c("max_posterior", "weighted")) {
        stop("posterior is required when mode = '", mode, "'. Use mode = 'density' to request the density-based fallback.", call. = FALSE)
    }
    if (!is.null(posterior) && mode %in% c("weighted", "max_posterior")) {
        posterior = prepare_grouped_posterior(
            posterior = posterior,
            n = n,
            target_cell_types = cell_types,
            cell_type_groups = cell_type_groups
        )
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

    if (mode %in% c("max_posterior", "density")) {
        if (total_only) {
            max_type = rep(cell_types, n)
        } else if (mode == "max_posterior") {
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

basis_watershed_transcript_cells <- function(
    basis_partition,
    transcripts_df,
    posterior = NULL,
    mode = c("max_posterior", "weighted", "density"),
    cell_type_groups = partition_cell_type_groups(basis_partition)
) {
    mode = match.arg(mode)
    if (is.null(basis_partition$transcript_basin) || is.null(basis_partition$basin_summary)) {
        stop("basis_partition must contain transcript_basin and basin_summary.", call. = FALSE)
    }
    cell_types = colnames(basis_partition$transcript_basin)
    n = nrow(basis_partition$transcript_basin)
    if (nrow(transcripts_df) != n) {
        stop("transcripts_df must have one row per transcript assignment.", call. = FALSE)
    }
    total_only = mode == "density" &&
        length(cell_types) == 1L &&
        isTRUE(basis_partition$parameters$mode == "total_density")

    if (is.null(posterior) && mode %in% c("max_posterior", "weighted")) {
        stop("posterior is required when mode = '", mode, "'. Use mode = 'density' to request the density-based fallback.", call. = FALSE)
    }
    if (!is.null(posterior)) {
        posterior = prepare_grouped_posterior(
            posterior = posterior,
            n = n,
            target_cell_types = cell_types,
            cell_type_groups = cell_type_groups
        )
    }

    basin_lookup = basin_cell_lookup(basis_partition)
    transcript_row = list()
    cell_out = list()
    weight_out = list()
    posterior_out = list()
    part_i = 0L

    append_rows = function(rows, cell_type, weight) {
        basin_id = basis_partition$transcript_basin[rows, cell_type]
        cell_name = basin_lookup[[cell_type]][as.character(basin_id)]
        keep = !is.na(cell_name)
        if (!any(keep)) {
            return(NULL)
        }
        rows = rows[keep]
        cell_name = cell_name[keep]
        weight = weight[keep]
        cell_type_posterior = rep(NA_real_, length(rows))
        if (!is.null(posterior) && cell_type %in% colnames(posterior)) {
            cell_type_posterior = posterior[rows, cell_type]
        }
        part_i <<- part_i + 1L
        transcript_row[[part_i]] <<- rows
        cell_out[[part_i]] <<- unname(cell_name)
        weight_out[[part_i]] <<- weight
        posterior_out[[part_i]] <<- cell_type_posterior
        NULL
    }

    if (mode == "weighted") {
        for (cell_type in cell_types) {
            if (!cell_type %in% colnames(posterior)) next
            rows = which(!is.na(basis_partition$transcript_basin[, cell_type]) & posterior[, cell_type] > 0)
            if (length(rows) == 0L) next
            append_rows(rows, cell_type, posterior[rows, cell_type])
        }
    } else {
        if (total_only) {
            max_type = rep(cell_types, n)
        } else if (mode == "max_posterior") {
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
            append_rows(rows, cell_type, rep(1, length(rows)))
        }
    }

    if (part_i == 0L) {
        return(data.frame(
            transcript_row = integer(),
            cell = character(),
            assignment_weight = numeric(),
            cell_type_posterior = numeric()
        ))
    }
    data.frame(
        transcript_row = unlist(transcript_row, use.names = FALSE),
        cell = unlist(cell_out, use.names = FALSE),
        assignment_weight = unlist(weight_out, use.names = FALSE),
        cell_type_posterior = unlist(posterior_out, use.names = FALSE),
        stringsAsFactors = FALSE
    )
}

estimate_basis_watershed_domain_measure <- function(basis_partition, density_fit) {
    na_out = data.frame(
        n_domain_simplexes_total = NA_integer_,
        domain_measure_total = NA_real_,
        domain_area_total = NA_real_,
        domain_volume_total = NA_real_,
        n_active_domain_simplexes = NA_integer_,
        active_domain_measure = NA_real_,
        active_domain_area = NA_real_,
        active_domain_volume = NA_real_
    )
    if (is.null(density_fit) || is.null(density_fit$domain) || is.null(density_fit$domain$basis_id) || is.null(density_fit$domain$volume)) {
        return(list(measure = NULL, domain_metadata = na_out))
    }
    if (is.null(density_fit$basis_points) || nrow(density_fit$basis_points) != nrow(basis_partition$basis_points)) {
        stop("density_fit must use the same density basis as basis_partition to estimate domain area/volume.", call. = FALSE)
    }
    domain_basis_id = as.matrix(density_fit$domain$basis_id)
    storage.mode(domain_basis_id) = "integer"
    n_simplex = nrow(domain_basis_id)
    simplex_volume = as.numeric(density_fit$domain$volume)
    if (length(simplex_volume) == 1L) {
        simplex_volume = rep(simplex_volume, n_simplex)
    }
    if (length(simplex_volume) != n_simplex) {
        stop("density_fit$domain$volume must have length 1 or one value per simplex.", call. = FALSE)
    }

    d = ncol(basis_partition$basis_points)
    raw = basis_watershed_domain_measure_cpp(
        domain_basis_id = domain_basis_id,
        simplex_volume = simplex_volume,
        active_start = basis_partition$active_start,
        type_density = basis_partition$type_density_basis,
        basin = basis_partition$basin
    )

    domain_metadata = raw$domain_metadata
    domain_metadata$domain_area_total = if (d == 2L) domain_metadata$domain_measure_total else NA_real_
    domain_metadata$domain_volume_total = if (d == 3L) domain_metadata$domain_measure_total else NA_real_
    domain_metadata$active_domain_area = if (d == 2L) domain_metadata$active_domain_measure else NA_real_
    domain_metadata$active_domain_volume = if (d == 3L) domain_metadata$active_domain_measure else NA_real_
    domain_metadata = domain_metadata[, names(na_out), drop = FALSE]

    measure = raw$measure
    if (nrow(measure) > 0L) {
        basin_key = paste(
            basis_partition$basin_summary$cell_type,
            basis_partition$basin_summary$basin,
            sep = "\r"
        )
        cell_type = basis_partition$active_cell_types[measure$type_index]
        cell_key = paste(cell_type, measure$basin, sep = "\r")
        cell = as.character(basis_partition$basin_summary$cell)[match(cell_key, basin_key)]
        keep = !is.na(cell)
        measure = data.frame(
            cell = cell[keep],
            n_domain_simplexes = measure$n_domain_simplexes[keep],
            basis_measure = measure$basis_measure[keep],
            stringsAsFactors = FALSE
        )
        measure$basis_area = if (d == 2L) measure$basis_measure else NA_real_
        measure$basis_volume = if (d == 3L) measure$basis_measure else NA_real_
    } else {
        measure = NULL
    }

    list(
        measure = measure,
        domain_metadata = domain_metadata
    )
}

estimate_basis_watershed_cell_measure <- function(basis_partition, density_fit) {
    estimate_basis_watershed_domain_measure(basis_partition, density_fit)$measure
}

summarize_basis_watershed_domain_measure <- function(basis_partition, density_fit) {
    estimate_basis_watershed_domain_measure(basis_partition, density_fit)$domain_metadata
}

#' Collect basis-watershed gene counts with cell and transcript metadata
#'
#' Returns the sparse gene-by-cell count matrix together with transcript-level
#' cell assignments and per-cell spatial metadata for plotting.
#'
#' @param basis_partition Result from [partition_basis_watershed_initial()] or
#'   [partition_basis_watershed_total()].
#' @param transcripts_df Transcript-level data frame corresponding to the
#'   density and segmentation fits.
#' @param posterior Optional transcript-by-cell-type posterior matrix. Required
#'   for `mode = "max_posterior"` and `mode = "weighted"`.
#' @param density_fit Optional density fit used to estimate basis-domain
#'   area/volume for each cell. Each active density-domain simplex is split
#'   fractionally across active cell types using the average active cell-type
#'   support of its vertices; each type-specific fraction is then assigned to
#'   that type's highest-density active basin for the simplex.
#' @param gene Character gene column name.
#' @param x,y,z Character coordinate column names.
#' @param mode Character count/assignment mode passed to
#'   [basis_watershed_gene_counts()].
#' @param cell_type_groups Optional named list mapping coarse watershed cell
#'   type names to fine posterior column names. Defaults to the mapping stored
#'   in `basis_partition`, if present.
#'
#' @return A list with `counts`, `cell_metadata`, `transcript_cells`, and
#'   `domain_metadata`. `domain_metadata` contains the total density-domain
#'   area/volume and the active-domain area/volume assigned across watershed
#'   cells.
#'   `transcript_cells` contains `transcript_row`, `cell`,
#'   `assignment_weight`, `cell_type_posterior`, the gene column, and available
#'   coordinate columns.
#' @export
basis_watershed_gene_count_data <- function(
    basis_partition,
    transcripts_df,
    posterior = NULL,
    density_fit = NULL,
    gene = "feature_name",
    x = "x_location",
    y = "y_location",
    z = "z_location",
    mode = c("max_posterior", "weighted", "density"),
    cell_type_groups = partition_cell_type_groups(basis_partition)
) {
    mode = match.arg(mode)
    counts = basis_watershed_gene_counts(
        basis_partition = basis_partition,
        transcripts_df = transcripts_df,
        posterior = posterior,
        gene = gene,
        mode = mode,
        cell_type_groups = cell_type_groups
    )
    transcript_cells = basis_watershed_transcript_cells(
        basis_partition = basis_partition,
        transcripts_df = transcripts_df,
        posterior = posterior,
        mode = mode,
        cell_type_groups = cell_type_groups
    )
    coord_cols = c(x, y, z)
    coord_cols = coord_cols[coord_cols %in% colnames(transcripts_df)]
    if (length(coord_cols) == 0L) {
        stop("No coordinate columns were found in transcripts_df.", call. = FALSE)
    }
    if (nrow(transcript_cells) > 0L) {
        transcript_cells = cbind(
            transcript_cells,
            transcripts_df[transcript_cells$transcript_row, c(gene, coord_cols), drop = FALSE]
        )
    } else {
        transcript_cells[[gene]] = character()
        for (coord_col in coord_cols) {
            transcript_cells[[coord_col]] = numeric()
        }
    }

    summary = basis_partition$basin_summary
    cell_metadata = summary[match(colnames(counts), as.character(summary$cell)), , drop = FALSE]
    rownames(cell_metadata) = NULL
    cell_metadata$cell = as.character(cell_metadata$cell)
    cell_metadata$n_detected_genes = Matrix::colSums(counts > 0)
    cell_metadata$count_sum = Matrix::colSums(counts)

    if (nrow(transcript_cells) > 0L) {
        assignment_n = rowsum(rep(1, nrow(transcript_cells)), transcript_cells$cell, reorder = FALSE)
        weight_sum = rowsum(transcript_cells$assignment_weight, transcript_cells$cell, reorder = FALSE)
        cell_metadata$n_transcript_assignments = as.integer(assignment_n[match(cell_metadata$cell, rownames(assignment_n)), 1L])
        cell_metadata$n_transcript_assignments[is.na(cell_metadata$n_transcript_assignments)] = 0L
        cell_metadata$assignment_weight_sum = as.numeric(weight_sum[match(cell_metadata$cell, rownames(weight_sum)), 1L])
        cell_metadata$assignment_weight_sum[is.na(cell_metadata$assignment_weight_sum)] = 0
        for (coord_col in coord_cols) {
            weighted_coord = rowsum(
                transcript_cells[[coord_col]] * transcript_cells$assignment_weight,
                transcript_cells$cell,
                reorder = FALSE
            )
            centroid = as.numeric(weighted_coord[match(cell_metadata$cell, rownames(weighted_coord)), 1L]) /
                cell_metadata$assignment_weight_sum
            centroid[!is.finite(centroid)] = NA_real_
            cell_metadata[[paste0(coord_col, "_centroid")]] = centroid
        }
    } else {
        cell_metadata$n_transcript_assignments = integer(nrow(cell_metadata))
        cell_metadata$assignment_weight_sum = numeric(nrow(cell_metadata))
        for (coord_col in coord_cols) {
            cell_metadata[[paste0(coord_col, "_centroid")]] = NA_real_
        }
    }

    domain_measure = estimate_basis_watershed_domain_measure(basis_partition, density_fit)
    measure = domain_measure$measure
    domain_metadata = domain_measure$domain_metadata
    if (!is.null(measure)) {
        cell_metadata = merge(cell_metadata, measure, by = "cell", all.x = TRUE, sort = FALSE)
        cell_metadata = cell_metadata[match(colnames(counts), cell_metadata$cell), , drop = FALSE]
        rownames(cell_metadata) = NULL
        no_measure = is.na(cell_metadata$basis_measure)
        cell_metadata$n_domain_simplexes[is.na(cell_metadata$n_domain_simplexes)] = 0L
        cell_metadata$basis_measure[no_measure] = 0
        if (ncol(basis_partition$basis_points) == 2L) {
            cell_metadata$basis_area[no_measure] = 0
        } else if (ncol(basis_partition$basis_points) == 3L) {
            cell_metadata$basis_volume[no_measure] = 0
        }
    } else {
        cell_metadata$n_domain_simplexes = NA_integer_
        cell_metadata$basis_measure = NA_real_
        cell_metadata$basis_area = NA_real_
        cell_metadata$basis_volume = NA_real_
    }

    list(
        counts = counts,
        cell_metadata = cell_metadata,
        transcript_cells = transcript_cells,
        domain_metadata = domain_metadata
    )
}
