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
    key = lo * (n + 1) + hi
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

normalize_target_cell_type_counts <- function(target_cell_type_counts, posterior) {
    if (is.null(target_cell_type_counts)) {
        return(NULL)
    }

    cell_types = colnames(posterior)
    if (is.null(cell_types)) {
        cell_types = as.character(seq_len(ncol(posterior)))
    }

    if (is.data.frame(target_cell_type_counts)) {
        if (!all(c("cell_type", "estimated_total_cells") %in% colnames(target_cell_type_counts))) {
            stop(
                "target_cell_type_counts data frame must contain cell_type and estimated_total_cells columns.",
                call. = FALSE
            )
        }
        target_df = target_cell_type_counts
        if ("estimable" %in% colnames(target_df)) {
            target_df = target_df[!is.na(target_df$estimable) & target_df$estimable, , drop = FALSE]
        }
        target = target_df$estimated_total_cells
        names(target) = as.character(target_df$cell_type)
    } else {
        target = as.numeric(target_cell_type_counts)
        if (is.null(names(target))) {
            if (length(target) != length(cell_types)) {
                stop(
                    "Unnamed target_cell_type_counts must have one value per posterior column.",
                    call. = FALSE
                )
            }
            names(target) = cell_types
        }
    }

    keep = names(target) %in% cell_types & is.finite(target) & target > 0
    target = target[keep]
    if (length(target) == 0L) {
        return(NULL)
    }
    target[cell_types[cell_types %in% names(target)]]
}

weighted_quantile_type_limit <- function(values, probs) {
    stats::quantile(values, probs = probs, names = FALSE, type = 7, na.rm = TRUE)
}

estimate_prior_cell_size_limits <- function(
    transcripts_df,
    posterior,
    cell_id = "cell_id",
    unassigned = "UNASSIGNED",
    min_transcripts_per_prior_cell = 10L,
    purity_threshold = 0.8,
    upper_quantile = 0.99,
    min_cells_per_type = 20L,
    shrink_to_global = TRUE
) {
    if (!cell_id %in% colnames(transcripts_df)) {
        stop("Column '", cell_id, "' not found in transcripts_df.", call. = FALSE)
    }
    posterior = normalize_posterior_matrix(posterior, nrow(transcripts_df))
    if (min_transcripts_per_prior_cell < 1L) {
        stop("min_transcripts_per_prior_cell must be at least 1.", call. = FALSE)
    }
    check_finite_scalar(purity_threshold, "purity_threshold", lower = 0)
    if (purity_threshold > 1) {
        stop("purity_threshold must be in [0, 1].", call. = FALSE)
    }
    check_finite_scalar(upper_quantile, "upper_quantile", lower = 0)
    if (upper_quantile <= 0 || upper_quantile >= 1) {
        stop("upper_quantile must be in (0, 1).", call. = FALSE)
    }
    if (min_cells_per_type < 1L) {
        stop("min_cells_per_type must be at least 1.", call. = FALSE)
    }

    assigned = !is.na(transcripts_df[[cell_id]]) & transcripts_df[[cell_id]] != unassigned
    if (!any(assigned)) {
        stop("No assigned prior cells were found.", call. = FALSE)
    }
    prior_counts = table(transcripts_df[[cell_id]][assigned])
    reliable_ids = names(prior_counts)[prior_counts >= min_transcripts_per_prior_cell]
    reliable = assigned & transcripts_df[[cell_id]] %in% reliable_ids
    if (!any(reliable)) {
        stop("No prior cells passed min_transcripts_per_prior_cell.", call. = FALSE)
    }

    cell_mass = rowsum(posterior[reliable, , drop = FALSE], group = transcripts_df[[cell_id]][reliable], reorder = FALSE)
    cell_n = as.numeric(table(transcripts_df[[cell_id]][reliable])[rownames(cell_mass)])
    cell_fraction = sweep(cell_mass, 1L, rowSums(cell_mass), "/")
    dominant = max.col(cell_fraction, ties.method = "first")
    purity = cell_fraction[cbind(seq_len(nrow(cell_fraction)), dominant)]
    cell_types = colnames(posterior)
    if (is.null(cell_types)) {
        cell_types = as.character(seq_len(ncol(posterior)))
    }

    keep = purity >= purity_threshold
    global_limit = weighted_quantile_type_limit(cell_n[keep], upper_quantile)
    if (!is.finite(global_limit)) {
        global_limit = weighted_quantile_type_limit(cell_n, upper_quantile)
    }

    limits = rep(global_limit, length(cell_types))
    n_cells = integer(length(cell_types))
    raw_limits = rep(NA_real_, length(cell_types))
    names(limits) = cell_types
    names(n_cells) = cell_types
    names(raw_limits) = cell_types

    for (k in seq_along(cell_types)) {
        idx = keep & dominant == k
        n_cells[k] = sum(idx)
        if (n_cells[k] > 0L) {
            raw_limits[k] = weighted_quantile_type_limit(cell_n[idx], upper_quantile)
            limits[k] = raw_limits[k]
        }
        if (isTRUE(shrink_to_global) && n_cells[k] < min_cells_per_type) {
            type_weight = n_cells[k] / min_cells_per_type
            limits[k] = type_weight * limits[k] + (1 - type_weight) * global_limit
        }
    }

    data.frame(
        cell_type = cell_types,
        max_transcripts = as.numeric(limits),
        raw_max_transcripts = as.numeric(raw_limits),
        high_purity_prior_cells = as.integer(n_cells),
        global_max_transcripts = as.numeric(global_limit),
        purity_threshold = purity_threshold,
        upper_quantile = upper_quantile,
        min_cells_per_type = as.integer(min_cells_per_type),
        shrink_to_global = isTRUE(shrink_to_global)
    )
}

normalize_cell_size_limits <- function(cell_size_limits, posterior) {
    if (is.null(cell_size_limits)) {
        return(NULL)
    }
    cell_types = colnames(posterior)
    if (is.null(cell_types)) {
        cell_types = as.character(seq_len(ncol(posterior)))
    }

    if (is.data.frame(cell_size_limits)) {
        if (!all(c("cell_type", "max_transcripts") %in% colnames(cell_size_limits))) {
            stop("cell_size_limits data frame must contain cell_type and max_transcripts columns.", call. = FALSE)
        }
        limits = cell_size_limits$max_transcripts
        names(limits) = as.character(cell_size_limits$cell_type)
    } else {
        limits = as.numeric(cell_size_limits)
        if (is.null(names(limits))) {
            if (length(limits) != length(cell_types)) {
                stop("Unnamed cell_size_limits must have one value per posterior column.", call. = FALSE)
            }
            names(limits) = cell_types
        }
    }

    keep = names(limits) %in% cell_types & is.finite(limits) & limits > 0
    limits = limits[keep]
    if (length(limits) == 0L) {
        return(NULL)
    }
    limits[cell_types[cell_types %in% names(limits)]]
}

target_count_loss <- function(current_counts, target_counts, target_count_sd) {
    shared = intersect(names(target_counts), names(current_counts))
    if (length(shared) == 0L) {
        return(NA_real_)
    }
    sum(((current_counts[shared] - target_counts[shared]) / target_count_sd[shared])^2)
}

posterior_weighted_component_counts <- function(
    component_mass,
    active,
    component_size = NULL,
    min_transcripts_per_counted_cell = 1L
) {
    counted = active
    if (!is.null(component_size)) {
        counted = counted & component_size >= min_transcripts_per_counted_cell
    }
    if (!any(counted)) {
        out = rep(0, ncol(component_mass))
        names(out) = colnames(component_mass)
        return(out)
    }
    active_mass = component_mass[counted, , drop = FALSE]
    component_total = rowSums(active_mass)
    component_fraction = sweep(active_mass, 1L, component_total, "/")
    colSums(component_fraction)
}

merge_watershed_basins <- function(
    basin,
    adjacency,
    saddle_ratio_threshold,
    posterior_js_threshold,
    min_transcripts_per_cell,
    basin_posterior_sum = NULL,
    target_cell_type_counts = NULL,
    target_count_sd = NULL,
    use_target_counts = TRUE,
    target_merge_saddle_ratio = 0.3,
    target_merge_posterior_js = 0.4,
    target_min_transcripts_per_cell = 1L,
    cell_size_limits = NULL,
    size_limit_slack = 1.25
) {
    if (nrow(adjacency) == 0L) {
        return(list(
            cell_assignment = basin,
            merge_history = data.frame(),
            count_trajectory = data.frame()
        ))
    }

    basin_levels = sort(unique(basin))
    basin_id = match(basin, basin_levels)
    uf = make_union_find(length(basin_levels))
    basin_size = tabulate(basin_id, nbins = length(basin_levels))
    active = rep(TRUE, length(basin_levels))

    target_counts = NULL
    current_counts = NULL
    current_count_loss = NA_real_
    count_trajectory = list()
    if (isTRUE(use_target_counts) && !is.null(target_cell_type_counts) && !is.null(basin_posterior_sum)) {
        target_counts = target_cell_type_counts
        if (is.null(target_count_sd)) {
            target_count_sd = sqrt(pmax(target_counts, 1))
        } else {
            target_count_sd = as.numeric(target_count_sd)
            if (is.null(names(target_count_sd))) {
                names(target_count_sd) = names(target_counts)
            }
            target_count_sd = target_count_sd[names(target_counts)]
        }
        if (any(!is.finite(target_count_sd)) || any(target_count_sd <= 0)) {
            stop("target_count_sd must contain positive finite values for every target cell type.", call. = FALSE)
        }

        current_counts = posterior_weighted_component_counts(
            basin_posterior_sum,
            active,
            basin_size,
            target_min_transcripts_per_cell
        )
        current_count_loss = target_count_loss(current_counts, target_counts, target_count_sd)
        count_trajectory[[1L]] = data.frame(
            merge_step = 0L,
            n_cells = sum(active),
            n_counted_cells = sum(active & basin_size >= target_min_transcripts_per_cell),
            count_loss = current_count_loss,
            as.data.frame(as.list(stats::setNames(current_counts, paste0("count_", names(current_counts)))))
        )
    }

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

        merged_size = basin_size[ra] + basin_size[rb]
        merged_mass = basin_posterior_sum[ra, , drop = FALSE] + basin_posterior_sum[rb, , drop = FALSE]
        size_guard_pass = TRUE
        merged_size_limit = NA_real_
        if (!is.null(cell_size_limits)) {
            shared_limits = intersect(names(cell_size_limits), colnames(basin_posterior_sum))
            if (length(shared_limits) > 0L) {
                merged_profile = merged_mass[, shared_limits, drop = TRUE] / sum(merged_mass)
                merged_size_limit = sum(merged_profile * cell_size_limits[shared_limits]) * size_limit_slack
                size_guard_pass = merged_size <= merged_size_limit
            }
        }

        weak_boundary = adjacency$saddle_ratio[r] >= saddle_ratio_threshold &&
            adjacency$basin_js[r] <= posterior_js_threshold
        small_basin = basin_size[ra] < min_transcripts_per_cell ||
            basin_size[rb] < min_transcripts_per_cell
        small_compatible = small_basin && adjacency$basin_js[r] <= posterior_js_threshold

        target_supported = FALSE
        count_loss_after = NA_real_
        if (!is.null(target_counts)) {
            plausible_for_target = adjacency$saddle_ratio[r] >= target_merge_saddle_ratio &&
                adjacency$basin_js[r] <= target_merge_posterior_js
            if (plausible_for_target) {
                shared = intersect(names(target_counts), colnames(basin_posterior_sum))
                counted_a = basin_size[ra] >= target_min_transcripts_per_cell
                counted_b = basin_size[rb] >= target_min_transcripts_per_cell
                counted_merged = (basin_size[ra] + basin_size[rb]) >= target_min_transcripts_per_cell
                counts_after = current_counts
                if (counted_a) {
                    before_a = basin_posterior_sum[ra, shared, drop = TRUE]
                    frac_a = before_a / sum(basin_posterior_sum[ra, ])
                    counts_after[shared] = counts_after[shared] - frac_a
                }
                if (counted_b) {
                    before_b = basin_posterior_sum[rb, shared, drop = TRUE]
                    frac_b = before_b / sum(basin_posterior_sum[rb, ])
                    counts_after[shared] = counts_after[shared] - frac_b
                }
                if (counted_merged) {
                    frac_merged = merged_mass[, shared, drop = TRUE] / sum(merged_mass)
                    counts_after[shared] = counts_after[shared] + frac_merged
                }
                count_loss_after = target_count_loss(counts_after, target_counts, target_count_sd)
                target_supported = is.finite(count_loss_after) && count_loss_after < current_count_loss
            }
        }

        should_merge = size_guard_pass && (weak_boundary || small_compatible || target_supported)

        if (should_merge) {
            reason = if (weak_boundary) {
                "weak_boundary"
            } else if (small_compatible) {
                "small_compatible_basin"
            } else {
                "target_count_improvement"
            }
            count_loss_before = current_count_loss
            new_root = uf$union(ra, rb)
            old_root = if (new_root == ra) rb else ra
            old_size = basin_size[ra] + basin_size[rb]
            basin_size[new_root] = old_size
            active[old_root] = FALSE

            if (!is.null(target_counts)) {
                basin_posterior_sum[new_root, ] = basin_posterior_sum[ra, ] + basin_posterior_sum[rb, ]
                if (new_root != ra) {
                    basin_posterior_sum[ra, ] = 0
                }
                if (new_root != rb) {
                    basin_posterior_sum[rb, ] = 0
                }
                current_counts = posterior_weighted_component_counts(
                    basin_posterior_sum,
                    active,
                    basin_size,
                    target_min_transcripts_per_cell
                )
                current_count_loss = target_count_loss(current_counts, target_counts, target_count_sd)
                count_trajectory[[length(count_trajectory) + 1L]] = data.frame(
                    merge_step = length(merge_rows) + 1L,
                    n_cells = sum(active),
                    n_counted_cells = sum(active & basin_size >= target_min_transcripts_per_cell),
                    count_loss = current_count_loss,
                    as.data.frame(as.list(stats::setNames(current_counts, paste0("count_", names(current_counts)))))
                )
            }

            merge_rows[[length(merge_rows) + 1L]] = data.frame(
                basin_a = adjacency$basin_a[r],
                basin_b = adjacency$basin_b[r],
                saddle_ratio = adjacency$saddle_ratio[r],
                basin_js = adjacency$basin_js[r],
                merged_transcripts = merged_size,
                merged_size_limit = merged_size_limit,
                count_loss_before = count_loss_before,
                count_loss_after = current_count_loss,
                reason = reason
            )
        }
    }

    roots = vapply(basin_id, uf$find, integer(1L))
    cell_assignment = match(roots, sort(unique(roots)))
    merge_history = if (length(merge_rows) == 0L) data.frame() else do.call(rbind, merge_rows)
    count_trajectory = if (length(count_trajectory) == 0L) data.frame() else do.call(rbind, count_trajectory)
    list(
        cell_assignment = cell_assignment,
        merge_history = merge_history,
        count_trajectory = count_trajectory
    )
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
#' @param use_target_counts Logical; if `TRUE`, use estimated or supplied
#'   cell-type target counts to accept additional locally plausible merges that
#'   reduce posterior-weighted count loss.
#' @param target_cell_type_counts Optional numeric vector or data frame of
#'   target cell counts by cell type. If `NULL`, the function uses the prior
#'   `cell_id` calibration when available.
#' @param target_count_sd Optional positive scale for each target count. If
#'   `NULL`, uses `sqrt(max(target, 1))`.
#' @param target_merge_saddle_ratio Minimum saddle ratio for target-count-guided
#'   merges.
#' @param target_merge_posterior_js Maximum basin posterior Jensen-Shannon
#'   divergence for target-count-guided merges.
#' @param target_min_transcripts_per_cell Minimum transcript count for watershed
#'   components to be counted in the target-count merge loss. If `NULL`, uses
#'   `prior_min_transcripts_per_cell` so prior target estimation and new-cell
#'   counting use the same transcript-count filter.
#' @param use_size_limits Logical; if `TRUE`, estimate or use type-specific
#'   upper transcript-count limits and block merges that would create
#'   implausibly large cells.
#' @param cell_size_limits Optional numeric vector or data frame with
#'   `cell_type` and `max_transcripts` columns. If `NULL`, limits are estimated
#'   from high-purity reliable prior cells when possible.
#' @param size_limit_purity_threshold Minimum dominant posterior fraction for a
#'   prior cell to contribute to size-limit estimation.
#' @param size_limit_upper_quantile Upper quantile of prior-cell transcript
#'   counts used as the per-type size limit.
#' @param size_limit_min_cells_per_type Minimum high-purity prior cells before a
#'   type-specific limit is trusted without shrinkage to the global limit.
#' @param size_limit_slack Multiplicative slack applied to type-weighted size
#'   limits during merging.
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
    use_target_counts = TRUE,
    target_cell_type_counts = NULL,
    target_count_sd = NULL,
    target_merge_saddle_ratio = 0.3,
    target_merge_posterior_js = 0.4,
    target_min_transcripts_per_cell = NULL,
    use_size_limits = TRUE,
    cell_size_limits = NULL,
    size_limit_purity_threshold = 0.8,
    size_limit_upper_quantile = 0.99,
    size_limit_min_cells_per_type = 20L,
    size_limit_slack = 1.25,
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
    check_finite_scalar(target_merge_saddle_ratio, "target_merge_saddle_ratio", lower = 0)
    check_finite_scalar(target_merge_posterior_js, "target_merge_posterior_js", lower = 0)
    check_finite_scalar(size_limit_purity_threshold, "size_limit_purity_threshold", lower = 0)
    if (size_limit_purity_threshold > 1) {
        stop("size_limit_purity_threshold must be in [0, 1].", call. = FALSE)
    }
    check_finite_scalar(size_limit_upper_quantile, "size_limit_upper_quantile", lower = 0)
    if (size_limit_upper_quantile <= 0 || size_limit_upper_quantile >= 1) {
        stop("size_limit_upper_quantile must be in (0, 1).", call. = FALSE)
    }
    check_finite_scalar(size_limit_slack, "size_limit_slack", lower = 0, lower_strict = TRUE)
    if (size_limit_min_cells_per_type < 1L) {
        stop("size_limit_min_cells_per_type must be at least 1.", call. = FALSE)
    }
    if (min_transcripts_per_cell < 1L) {
        stop("min_transcripts_per_cell must be at least 1.", call. = FALSE)
    }
    if (is.null(target_min_transcripts_per_cell)) {
        target_min_transcripts_per_cell = prior_min_transcripts_per_cell
    }
    if (
        length(target_min_transcripts_per_cell) != 1L ||
        !is.finite(target_min_transcripts_per_cell) ||
        target_min_transcripts_per_cell < 1 ||
        target_min_transcripts_per_cell != as.integer(target_min_transcripts_per_cell)
    ) {
        stop("target_min_transcripts_per_cell must be NULL or a positive integer.", call. = FALSE)
    }
    target_min_transcripts_per_cell = as.integer(target_min_transcripts_per_cell)

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
    if (is.null(target_cell_type_counts)) {
        target_cell_type_counts = prior_counts
    }
    target_counts = normalize_target_cell_type_counts(target_cell_type_counts, posterior)
    prior_size_limits = NULL
    if (isTRUE(use_size_limits) && is.null(cell_size_limits) && cell_id %in% colnames(transcripts_df)) {
        prior_size_limits = tryCatch(
            estimate_prior_cell_size_limits(
                transcripts_df = transcripts_df,
                posterior = posterior,
                cell_id = cell_id,
                unassigned = unassigned,
                min_transcripts_per_prior_cell = prior_min_transcripts_per_cell,
                purity_threshold = size_limit_purity_threshold,
                upper_quantile = size_limit_upper_quantile,
                min_cells_per_type = size_limit_min_cells_per_type,
                shrink_to_global = TRUE
            ),
            error = function(e) {
                if (show_progress) {
                    message("Skipping prior cell size-limit calibration: ", conditionMessage(e))
                }
                NULL
            }
        )
    }
    if (is.null(cell_size_limits)) {
        cell_size_limits = prior_size_limits
    }
    normalized_size_limits = if (isTRUE(use_size_limits)) {
        normalize_cell_size_limits(cell_size_limits, posterior)
    } else {
        NULL
    }
    basin_posterior_sum = rowsum(posterior, group = initial_basin, reorder = TRUE)

    merge = merge_watershed_basins(
        basin = initial_basin,
        adjacency = basin_adjacency,
        saddle_ratio_threshold = saddle_ratio_threshold,
        posterior_js_threshold = posterior_js_threshold,
        min_transcripts_per_cell = min_transcripts_per_cell,
        basin_posterior_sum = basin_posterior_sum,
        target_cell_type_counts = target_counts,
        target_count_sd = target_count_sd,
        use_target_counts = use_target_counts,
        target_merge_saddle_ratio = target_merge_saddle_ratio,
        target_merge_posterior_js = target_merge_posterior_js,
        target_min_transcripts_per_cell = target_min_transcripts_per_cell,
        cell_size_limits = normalized_size_limits,
        size_limit_slack = size_limit_slack
    )

    final_cell = merge$cell_assignment
    df = transcripts_df
    df$watershed_basin = initial_basin
    df$watershed_cell = final_cell
    df$watershed_density = density

    list(
        transcripts_df = df,
        cell_assignment = final_cell,
        initial_basin = initial_basin,
        parent = parent,
        density = density,
        edges = transform(edges, posterior_js = edge_js),
        basin_adjacency = basin_adjacency,
        merge_history = merge$merge_history,
        count_trajectory = merge$count_trajectory,
        target_cell_type_counts = target_counts,
        cell_size_limits = normalized_size_limits,
        prior_cell_size_limits = prior_size_limits,
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
            min_transcripts_per_cell = min_transcripts_per_cell,
            use_target_counts = use_target_counts,
            target_merge_saddle_ratio = target_merge_saddle_ratio,
            target_merge_posterior_js = target_merge_posterior_js,
            target_min_transcripts_per_cell = target_min_transcripts_per_cell,
            use_size_limits = use_size_limits,
            size_limit_purity_threshold = size_limit_purity_threshold,
            size_limit_upper_quantile = size_limit_upper_quantile,
            size_limit_min_cells_per_type = size_limit_min_cells_per_type,
            size_limit_slack = size_limit_slack
        )
    )
}
