#' Profile Memory Usage Across `run_crf()` Stages
#'
#' Profiles the major stages of [run_crf()] on a given input and reports
#' stage-level memory usage, output dimensions, and object sizes.
#'
#' When the optional `peakRAM` package is installed, this function records peak
#' and total RAM usage for each stage. Otherwise it falls back to reporting
#' stage output sizes and allocation estimates based on `Rprofmem()`.
#'
#' @param transcripts_df A data frame containing transcript-level data.
#' @param cell_signatures A numeric matrix (gene x cell types) of probabilities
#'   or likelihoods for each gene within a cell type.
#' @param x Character; column name for x-coordinates.
#' @param y Character; column name for y-coordinates.
#' @param z Character; column name for z-coordinates.
#' @param gene Character; column name for gene identifiers.
#' @param qv Character; column name for quality values.
#' @param is_gene Character; column name marking rows to retain as genes.
#' @param qv_threshold Numeric; minimum quality value to retain transcripts.
#' @param n_neighbors Integer; number of nearest neighbors to connect.
#' @param dist_threshold Numeric; maximum spatial distance for edge inclusion.
#' @param same_label_ratio Numeric; same-label ratio used for Potts edge weights.
#' @param max_iter Integer; maximum LBP iterations passed to
#'   [potts_lbp_parallel_cpp()].
#' @param damping Numeric; damping factor passed to [potts_lbp_parallel_cpp()].
#' @param tol Numeric; convergence tolerance passed to [potts_lbp_parallel_cpp()].
#' @param num_cores Integer; number of CPU cores to use for parallel inference.
#' @param capture_objects Logical; if `TRUE`, include stage outputs in the
#'   returned object. Set to `FALSE` to avoid retaining large intermediates.
#' @param use_peakRAM Logical; if `TRUE` and the `peakRAM` package is installed,
#'   record stage peak memory usage.
#'
#' @return A list with components:
#'   \item{summary}{Data frame with one row per stage and columns describing
#'     memory usage and output dimensions.}
#'   \item{lbp_estimate_mb}{Approximate lower-bound memory footprint for the
#'     C++ LBP working buffers and returned marginals.}
#'   \item{used_peakRAM}{Logical; whether `peakRAM` was used.}
#'   \item{stages}{Named list of per-stage details and, optionally, outputs.}
#'
#' @examples
#' \dontrun{
#' prof <- profile_run_crf(transcripts_df, cell_signatures, capture_objects = FALSE)
#' prof$summary
#' prof$lbp_estimate_mb
#' }
#'
#' @export
profile_run_crf <- function(
    transcripts_df,
    cell_signatures,
    x = "x_location",
    y = "y_location",
    z = "z_location",
    gene = "feature_name",
    qv = "qv",
    is_gene = "is_gene",
    qv_threshold = 20,
    n_neighbors = 10,
    dist_threshold = 2,
    same_label_ratio = 5,
    max_iter = 50L,
    damping = 1.0,
    tol = 0.1,
    num_cores = 6,
    capture_objects = FALSE,
    use_peakRAM = TRUE
) {

    use_peakRAM = isTRUE(use_peakRAM) && requireNamespace("peakRAM", quietly = TRUE)
    strip_stage_value = function(entry) {
        entry$value = NULL
        entry
    }

    run_stage = function(stage_name, expr, info = NULL) {
        gc(verbose = FALSE)
        expr_sub = substitute(expr)
        eval_env = parent.frame()
        value = NULL
        evaluate_stage = function() {
            value <<- eval(expr_sub, envir = eval_env)
            value
        }

        if (use_peakRAM) {
            peak = peakRAM::peakRAM(evaluate_stage())
            total_alloc_mb = peak[["Total_RAM_Used_MiB"]][1]
            peak_alloc_mb = peak[["Peak_RAM_Used_MiB"]][1]
            alloc_source = "peakRAM"
        } else {
            tmp = tempfile("profile_run_crf_", fileext = ".out")
            utils::Rprofmem(tmp)
            value = evaluate_stage()
            utils::Rprofmem(NULL)
            mem_lines = readLines(tmp, warn = FALSE)
            unlink(tmp)

            alloc_bytes = suppressWarnings(sum(as.numeric(sub(" .*", "", mem_lines)), na.rm = TRUE))
            if (!is.finite(alloc_bytes)) {
                alloc_bytes = NA_real_
            }

            total_alloc_mb = alloc_bytes / (1024 ^ 2)
            peak_alloc_mb = NA_real_
            alloc_source = "Rprofmem"
        }

        size_mb = as.numeric(utils::object.size(value)) / (1024 ^ 2)
        entry = list(
            stage = stage_name,
            value = value,
            size_mb = size_mb,
            total_alloc_mb = total_alloc_mb,
            peak_alloc_mb = peak_alloc_mb,
            alloc_source = alloc_source,
            info = info
        )

        entry
    }

    stage_filter = run_stage("filter", {
        df = transcripts_df

        if (qv %in% colnames(df)) {
            df = df[df[[qv]] >= qv_threshold, , drop = FALSE]
        }

        if (is_gene %in% colnames(df)) {
            keep_gene = df[[is_gene]]
            keep_gene[is.na(keep_gene)] = FALSE
            df = df[keep_gene, , drop = FALSE]
        }

        keep_signature = df[[gene]] %in% rownames(cell_signatures)
        df[keep_signature, , drop = FALSE]
    })
    filtered_df = stage_filter$value

    stage_hilbert = run_stage("hilbert_sort", {
        idx = hilbert_index(as.matrix(filtered_df[, c(x, y), drop = FALSE]))
        filtered_df[order(idx), , drop = FALSE]
    }, info = list(n_nodes = nrow(filtered_df)))
    ordered_df = stage_hilbert$value

    stage_knn = run_stage("knn", {
        RANN::nn2(ordered_df[, c(x, y, z), drop = FALSE], k = n_neighbors + 1)
    }, info = list(n_nodes = nrow(ordered_df), k = n_neighbors + 1))
    knn = stage_knn$value

    edge_mask = as.numeric(knn$nn.dists[, -1, drop = FALSE]) <= dist_threshold
    n_directed_edges = sum(edge_mask)

    stage_graph = run_stage("graph_build", {
        from = rep(seq_len(nrow(knn$nn.idx)), n_neighbors)
        to = as.numeric(knn$nn.idx[, -1, drop = FALSE])
        keep = as.numeric(knn$nn.dists[, -1, drop = FALSE]) <= dist_threshold

        adj = Matrix::sparseMatrix(
            i = c(from[keep], to[keep]),
            j = c(to[keep], from[keep]),
            x = 1,
            dims = c(nrow(ordered_df), nrow(ordered_df))
        )
        adj@x = rep(log(same_label_ratio), length(adj@x))

        build_potts_lbp_graph(adj)
    }, info = list(n_nodes = nrow(ordered_df), n_directed_edges = 2L * n_directed_edges))
    graph = stage_graph$value

    stage_node_potentials = run_stage("node_potentials", {
        log(cell_signatures[ordered_df[[gene]], , drop = FALSE])
    }, info = list(n_nodes = nrow(ordered_df), n_labels = ncol(cell_signatures)))
    node_potentials = stage_node_potentials$value

    stage_lbp = run_stage("lbp", {
        potts_lbp_parallel_cpp(
            graph$adj_ptr,
            graph$adj_idx,
            graph$rev_idx,
            graph$edge_weights,
            node_potentials,
            max_iter = as.integer(max_iter),
            damping = damping,
            tol = tol,
            n_threads = min(num_cores, parallel::detectCores())
        )
    }, info = list(
        n_nodes = nrow(node_potentials),
        n_labels = ncol(node_potentials),
        n_directed_edges = length(graph$adj_idx)
    ))

    summary = do.call(
        rbind,
        lapply(
            list(stage_filter, stage_hilbert, stage_knn, stage_graph, stage_node_potentials, stage_lbp),
            function(entry) {
                info = entry$info
                data.frame(
                    stage = entry$stage,
                    object_size_mb = entry$size_mb,
                    total_alloc_mb = entry$total_alloc_mb,
                    peak_alloc_mb = entry$peak_alloc_mb,
                    alloc_source = entry$alloc_source,
                    n_nodes = if (is.null(info$n_nodes)) NA_integer_ else info$n_nodes,
                    n_labels = if (is.null(info$n_labels)) NA_integer_ else info$n_labels,
                    n_directed_edges = if (is.null(info$n_directed_edges)) NA_integer_ else info$n_directed_edges,
                    stringsAsFactors = FALSE
                )
            }
        )
    )

    lbp_estimate_mb = estimate_lbp_memory_mb(
        N = nrow(node_potentials),
        K = ncol(node_potentials),
        E = length(graph$adj_idx)
    )

    list(
        summary = summary,
        lbp_estimate_mb = lbp_estimate_mb,
        used_peakRAM = use_peakRAM,
        stages = list(
            filter = if (isTRUE(capture_objects)) stage_filter else strip_stage_value(stage_filter),
            hilbert_sort = if (isTRUE(capture_objects)) stage_hilbert else strip_stage_value(stage_hilbert),
            knn = if (isTRUE(capture_objects)) stage_knn else strip_stage_value(stage_knn),
            graph_build = if (isTRUE(capture_objects)) stage_graph else strip_stage_value(stage_graph),
            node_potentials = if (isTRUE(capture_objects)) stage_node_potentials else strip_stage_value(stage_node_potentials),
            lbp = if (isTRUE(capture_objects)) stage_lbp else strip_stage_value(stage_lbp)
        )
    )
}

estimate_lbp_memory_mb <- function(N, K, E) {
    bytes = E * K * 4 + E * K * 4 + N * K * 4 + N * K * 8
    bytes / (1024 ^ 2)
}
