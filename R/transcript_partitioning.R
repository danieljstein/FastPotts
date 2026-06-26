check_finite_scalar <- function(x, name, lower = -Inf, lower_strict = FALSE) {
    ok = length(x) == 1L && is.finite(x)
    if (ok) {
        ok = if (lower_strict) x > lower else x >= lower
    }
    if (!ok) {
        stop(name, " must be a finite scalar", if (is.finite(lower)) paste0(" >= ", lower) else "", ".", call. = FALSE)
    }
    invisible(NULL)
}

normalize_posterior_matrix <- function(posterior, n) {
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
    row_total = rowSums(posterior)
    if (any(row_total <= 0)) {
        stop("Every posterior row must have positive mass.", call. = FALSE)
    }
    sweep(posterior, 1L, row_total, "/")
}

posterior_js_divergence_edges <- function(posterior, from, to, eps = 1e-12) {
    posterior_js_divergence_edges_cpp(posterior, from, to, eps)
}

build_transcript_knn_edges <- function(coords, n_neighbors, max_distance = Inf) {
    n = nrow(coords)
    if (n < 2L) {
        stop("At least two transcripts are required to build a neighborhood graph.", call. = FALSE)
    }
    n_neighbors = min(as.integer(n_neighbors), n - 1L)
    nn = RANN::nn2(coords, k = n_neighbors + 1L)

    from = rep(seq_len(n), times = n_neighbors)
    to = as.vector(nn$nn.idx[, -1L, drop = FALSE])
    distance = as.vector(nn$nn.dists[, -1L, drop = FALSE])

    keep = to > 0L & from != to & distance <= max_distance
    from = from[keep]
    to = to[keep]
    distance = distance[keep]

    lo = pmin(from, to)
    hi = pmax(from, to)
    key = paste(lo, hi, sep = ":")
    ord = order(key, distance)
    lo = lo[ord]
    hi = hi[ord]
    distance = distance[ord]
    key = key[ord]
    keep_first = !duplicated(key)

    data.frame(
        from = lo[keep_first],
        to = hi[keep_first],
        distance = distance[keep_first]
    )
}

estimate_graph_density <- function(edges, posterior, bandwidth, mode = c("type_weighted", "total")) {
    mode = match.arg(mode)
    mode_id = match(mode, c("type_weighted", "total")) - 1L
    estimate_graph_density_cpp(edges$from, edges$to, edges$distance, posterior, bandwidth, mode_id)
}

find_density_ascent_parents <- function(edges, density, posterior_js, distance_weight, posterior_weight) {
    density_ascent_partition_cpp(
        edges$from,
        edges$to,
        edges$distance,
        density,
        posterior_js,
        distance_weight,
        posterior_weight
    )$parent
}

compress_density_roots <- function(parent) {
    root = parent
    for (i in seq_along(root)) {
        path = integer()
        j = i
        while (root[j] != j) {
            path = c(path, j)
            j = root[j]
        }
        if (length(path) > 0L) {
            root[path] = j
        }
        root[i] = j
    }
    match(root, sort(unique(root)))
}

summarize_cells <- function(cell_assignment, density, posterior) {
    cells = sort(unique(cell_assignment))
    K = ncol(posterior)
    posterior_sum = rowsum(posterior, group = cell_assignment, reorder = TRUE)
    posterior_mean = sweep(posterior_sum, 1L, rowSums(posterior_sum), "/")
    cell_type = max.col(posterior_sum, ties.method = "first")
    out = data.frame(
        cell_index = cells,
        n_transcripts = as.integer(tabulate(match(cell_assignment, cells), nbins = length(cells))),
        mode_density = as.numeric(tapply(density, cell_assignment, max)[as.character(cells)]),
        cell_type = if (is.null(colnames(posterior))) cell_type else colnames(posterior)[cell_type],
        posterior_mass = rowSums(posterior_sum)
    )
    type_cols = paste0("posterior_", if (is.null(colnames(posterior))) seq_len(K) else colnames(posterior))
    posterior_mean_df = as.data.frame(posterior_mean)
    colnames(posterior_mean_df) = type_cols
    cbind(out, posterior_mean_df)
}

build_basin_adjacency <- function(edges, basin, density, posterior, posterior_js) {
    build_basin_adjacency_cpp(
        edges$from,
        edges$to,
        edges$distance,
        basin,
        density,
        posterior,
        posterior_js
    )
}

make_union_find <- function(n) {
    parent = seq_len(n)
    size = rep(1L, n)

    find = function(x) {
        while (parent[x] != x) {
            parent[x] <<- parent[parent[x]]
            x = parent[x]
        }
        x
    }

    union = function(a, b) {
        ra = find(a)
        rb = find(b)
        if (ra == rb) {
            return(ra)
        }
        if (size[ra] < size[rb]) {
            tmp = ra
            ra = rb
            rb = tmp
        }
        parent[rb] <<- ra
        size[ra] <<- size[ra] + size[rb]
        ra
    }

    list(find = find, union = union)
}

merge_watershed_basins <- function(basin, adjacency, saddle_ratio_threshold, posterior_js_threshold, min_transcripts_per_cell) {
    if (nrow(adjacency) == 0L) {
        return(list(cell_assignment = basin, merge_history = data.frame()))
    }

    basin_levels = sort(unique(basin))
    basin_id = match(basin, basin_levels)
    uf = make_union_find(length(basin_levels))
    basin_size = tabulate(basin_id, nbins = length(basin_levels))

    adjacency = adjacency[order(-adjacency$saddle_ratio, adjacency$basin_js), , drop = FALSE]
    merge_rows = list()

    for (r in seq_len(nrow(adjacency))) {
        a = match(adjacency$basin_a[r], basin_levels)
        b = match(adjacency$basin_b[r], basin_levels)
        ra = uf$find(a)
        rb = uf$find(b)
        if (ra == rb) {
            next
        }

        weak_boundary = adjacency$saddle_ratio[r] >= saddle_ratio_threshold &&
            adjacency$basin_js[r] <= posterior_js_threshold
        small_basin = basin_size[ra] < min_transcripts_per_cell ||
            basin_size[rb] < min_transcripts_per_cell
        should_merge = weak_boundary || (small_basin && adjacency$basin_js[r] <= posterior_js_threshold)

        if (should_merge) {
            new_root = uf$union(ra, rb)
            old_size = basin_size[ra] + basin_size[rb]
            basin_size[new_root] = old_size
            merge_rows[[length(merge_rows) + 1L]] = data.frame(
                basin_a = adjacency$basin_a[r],
                basin_b = adjacency$basin_b[r],
                saddle_ratio = adjacency$saddle_ratio[r],
                basin_js = adjacency$basin_js[r],
                reason = if (weak_boundary) "weak_boundary" else "small_compatible_basin"
            )
        }
    }

    roots = vapply(basin_id, uf$find, integer(1L))
    cell_assignment = match(roots, sort(unique(roots)))
    merge_history = if (length(merge_rows) == 0L) data.frame() else do.call(rbind, merge_rows)
    list(cell_assignment = cell_assignment, merge_history = merge_history)
}

#' Estimate cell-type-specific cell counts from prior cell IDs
#'
#' Uses an existing transcript-level `cell_id` column as a weak calibration
#' target for how many cells of each posterior type may be present after adding
#' unassigned transcripts.
#'
#' @param transcripts_df Transcript-level data frame.
#' @param posterior Numeric matrix with one row per transcript and one column
#'   per cell type.
#' @param cell_id Character; prior cell ID column.
#' @param unassigned Character value used for transcripts without an assigned
#'   prior cell.
#' @param min_transcripts_per_prior_cell Integer minimum transcript count for a
#'   prior cell to be treated as reliable.
#'
#' @return A data frame with posterior-weighted assigned cell counts,
#'   transcript masses, and extrapolated total cell counts by type.
#' @export
estimate_prior_cell_type_counts <- function(
    transcripts_df,
    posterior,
    cell_id = "cell_id",
    unassigned = "UNASSIGNED",
    min_transcripts_per_prior_cell = 10L
) {
    if (!cell_id %in% colnames(transcripts_df)) {
        stop("Column '", cell_id, "' not found in transcripts_df.", call. = FALSE)
    }
    posterior = normalize_posterior_matrix(posterior, nrow(transcripts_df))
    if (min_transcripts_per_prior_cell < 1L) {
        stop("min_transcripts_per_prior_cell must be at least 1.", call. = FALSE)
    }

    assigned = !is.na(transcripts_df[[cell_id]]) & transcripts_df[[cell_id]] != unassigned
    prior_ids = transcripts_df[[cell_id]][assigned]
    if (!any(assigned)) {
        stop("No assigned prior cells were found.", call. = FALSE)
    }

    prior_counts = table(prior_ids)
    reliable_ids = names(prior_counts)[prior_counts >= min_transcripts_per_prior_cell]
    reliable = assigned & transcripts_df[[cell_id]] %in% reliable_ids
    if (!any(reliable)) {
        stop("No prior cells passed min_transcripts_per_prior_cell.", call. = FALSE)
    }

    cell_mass = rowsum(posterior[reliable, , drop = FALSE], group = transcripts_df[[cell_id]][reliable], reorder = FALSE)
    cell_type_fraction = sweep(cell_mass, 1L, rowSums(cell_mass), "/")

    total_mass = colSums(posterior)
    assigned_mass = colSums(cell_mass)
    assigned_cells = colSums(cell_type_fraction)
    estimable = assigned_mass > 0
    estimated_total = rep(NA_real_, length(assigned_mass))
    estimated_total[estimable] = assigned_cells[estimable] *
        total_mass[estimable] / assigned_mass[estimable]

    data.frame(
        cell_type = if (is.null(colnames(posterior))) seq_len(ncol(posterior)) else colnames(posterior),
        assigned_cells = as.numeric(assigned_cells),
        total_transcript_mass = as.numeric(total_mass),
        assigned_transcript_mass = as.numeric(assigned_mass),
        estimated_total_cells = as.numeric(estimated_total),
        estimable = estimable,
        reliable_prior_cells = nrow(cell_mass),
        min_transcripts_per_prior_cell = as.integer(min_transcripts_per_prior_cell)
    )
}

#' Partition transcripts into putative cells by graph watershed
#'
#' Builds a transcript neighborhood graph, estimates local transcript density,
#' assigns transcripts to density modes by graph ascent, and optionally merges
#' weakly separated basins using density-saddle and posterior-compatibility
#' diagnostics.
#'
#' @param transcripts_df Transcript-level data frame.
#' @param posterior Numeric matrix of transcript-by-cell-type posterior
#'   probabilities, such as `fit$marginals` from
#'   [spatial_basis_segmentation()].
#' @param x,y,z Character coordinate column names.
#' @param use_z Logical; if `TRUE`, include `z` in the spatial graph.
#' @param n_neighbors Integer number of nearest neighbors for the transcript
#'   graph.
#' @param max_distance Numeric maximum graph edge distance.
#' @param density_bandwidth Numeric kernel bandwidth. If `NULL`, uses the
#'   median graph edge distance.
#' @param density_mode Character; `"type_weighted"` uses posterior-weighted
#'   type density, while `"total"` uses total transcript density.
#' @param distance_weight Non-negative penalty for long uphill ascent edges.
#' @param posterior_weight Non-negative penalty for posterior divergence during
#'   uphill ascent.
#' @param saddle_ratio_threshold Merge adjacent basins when their connecting
#'   saddle is at least this fraction of the smaller mode density.
#' @param posterior_js_threshold Maximum Jensen-Shannon divergence allowed for
#'   compatible basin merges.
#' @param min_transcripts_per_cell Basins smaller than this are merged into a
#'   compatible neighbor when possible.
#' @param cell_id Character prior cell ID column used only for count
#'   calibration diagnostics.
#' @param unassigned Character prior cell ID value for unassigned transcripts.
#' @param prior_min_transcripts_per_cell Minimum transcript count for reliable
#'   prior cells in calibration diagnostics.
#' @param show_progress Logical; print progress messages.
#'
#' @return A list containing transcript assignments, density values, graph
#'   edges, initial basins, basin adjacency diagnostics, merge history, final
#'   cell summaries, and optional prior cell-type count calibration.
#' @export
partition_transcripts_watershed <- function(
    transcripts_df,
    posterior,
    x = "x_location",
    y = "y_location",
    z = "z_location",
    use_z = z %in% colnames(transcripts_df),
    n_neighbors = 20L,
    max_distance = Inf,
    density_bandwidth = NULL,
    density_mode = c("type_weighted", "total"),
    distance_weight = 0,
    posterior_weight = 1,
    saddle_ratio_threshold = 0.7,
    posterior_js_threshold = 0.2,
    min_transcripts_per_cell = 5L,
    cell_id = "cell_id",
    unassigned = "UNASSIGNED",
    prior_min_transcripts_per_cell = 10L,
    show_progress = TRUE
) {
    density_mode = match.arg(density_mode)
    n = nrow(transcripts_df)
    posterior = normalize_posterior_matrix(posterior, n)

    coord_cols = if (isTRUE(use_z)) c(x, y, z) else c(x, y)
    missing_cols = setdiff(coord_cols, colnames(transcripts_df))
    if (length(missing_cols) > 0L) {
        stop("Missing coordinate column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    coords = as.matrix(transcripts_df[, coord_cols, drop = FALSE])
    storage.mode(coords) = "double"
    if (any(!is.finite(coords))) {
        stop("Coordinate columns must contain finite values.", call. = FALSE)
    }

    if (n_neighbors < 1L) {
        stop("n_neighbors must be at least 1.", call. = FALSE)
    }
    if (length(max_distance) != 1L || is.na(max_distance) || max_distance <= 0) {
        stop("max_distance must be a positive numeric scalar or Inf.", call. = FALSE)
    }
    check_finite_scalar(distance_weight, "distance_weight", lower = 0)
    check_finite_scalar(posterior_weight, "posterior_weight", lower = 0)
    check_finite_scalar(saddle_ratio_threshold, "saddle_ratio_threshold", lower = 0)
    check_finite_scalar(posterior_js_threshold, "posterior_js_threshold", lower = 0)
    if (min_transcripts_per_cell < 1L) {
        stop("min_transcripts_per_cell must be at least 1.", call. = FALSE)
    }

    if (show_progress) {
        message("Building transcript neighborhood graph...")
    }
    edges = build_transcript_knn_edges(coords, n_neighbors = n_neighbors, max_distance = max_distance)
    if (nrow(edges) == 0L) {
        stop("No transcript graph edges were retained.", call. = FALSE)
    }

    if (is.null(density_bandwidth)) {
        density_bandwidth = stats::median(edges$distance)
    }
    check_finite_scalar(density_bandwidth, "density_bandwidth", lower = 0, lower_strict = TRUE)

    if (show_progress) {
        message("Estimating local transcript density...")
    }
    density = estimate_graph_density(edges, posterior, bandwidth = density_bandwidth, mode = density_mode)
    edge_js = posterior_js_divergence_edges(posterior, edges$from, edges$to)

    if (show_progress) {
        message("Assigning transcripts to density modes...")
    }
    ascent = density_ascent_partition_cpp(
        edges$from,
        edges$to,
        edges$distance,
        density,
        edge_js,
        distance_weight,
        posterior_weight
    )
    parent = ascent$parent
    initial_basin = ascent$basin

    if (show_progress) {
        message("Computing basin saddles and merge diagnostics...")
    }
    basin_adjacency = build_basin_adjacency(edges, initial_basin, density, posterior, edge_js)
    merge = merge_watershed_basins(
        basin = initial_basin,
        adjacency = basin_adjacency,
        saddle_ratio_threshold = saddle_ratio_threshold,
        posterior_js_threshold = posterior_js_threshold,
        min_transcripts_per_cell = min_transcripts_per_cell
    )

    final_cell = merge$cell_assignment
    df = transcripts_df
    df$watershed_basin = initial_basin
    df$watershed_cell = final_cell
    df$watershed_density = density

    prior_counts = NULL
    if (cell_id %in% colnames(transcripts_df)) {
        prior_counts = tryCatch(
            estimate_prior_cell_type_counts(
                transcripts_df = transcripts_df,
                posterior = posterior,
                cell_id = cell_id,
                unassigned = unassigned,
                min_transcripts_per_prior_cell = prior_min_transcripts_per_cell
            ),
            error = function(e) {
                if (show_progress) {
                    message("Skipping prior cell count calibration: ", conditionMessage(e))
                }
                NULL
            }
        )
    }

    list(
        transcripts_df = df,
        cell_assignment = final_cell,
        initial_basin = initial_basin,
        parent = parent,
        density = density,
        edges = transform(edges, posterior_js = edge_js),
        basin_adjacency = basin_adjacency,
        merge_history = merge$merge_history,
        cells = summarize_cells(final_cell, density, posterior),
        initial_basins = summarize_cells(initial_basin, density, posterior),
        prior_cell_type_counts = prior_counts,
        parameters = list(
            n_neighbors = n_neighbors,
            max_distance = max_distance,
            density_bandwidth = density_bandwidth,
            density_mode = density_mode,
            distance_weight = distance_weight,
            posterior_weight = posterior_weight,
            saddle_ratio_threshold = saddle_ratio_threshold,
            posterior_js_threshold = posterior_js_threshold,
            min_transcripts_per_cell = min_transcripts_per_cell
        )
    )
}
