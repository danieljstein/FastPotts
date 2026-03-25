#' Cell Type Assignment via Conditional Random Field
#'
#' Assigns cell type labels to individual transcripts using a Potts model
#' conditional random field (CRF) with loopy belief propagation inference.
#'
#' This function performs spatial cell type assignment by:
#' 1. Filtering transcripts by quality threshold
#' 2. Computing Hilbert-curve ordering for spatial coherence
#' 3. Building a k-nearest neighbor graph with spatial distance constraints
#' 4. Running loopy belief propagation on the Potts model to infer cell types
#'
#' The graph structure encodes both spatial proximity and the prior probability
#' of transcripts sharing the same cell type label.
#'
#' @param transcripts_df A data frame containing transcript-level data with
#'   columns for coordinates, gene names, and quality values.
#' @param cell_signatures A numeric matrix (gene x cell types) of probabilities
#'   or likelihoods for each gene within a cell type. Columns should sum to 1.
#' @param x Character; column name for x-coordinates (default: "x_location").
#' @param y Character; column name for y-coordinates (default: "y_location").
#' @param z Character; column name for z-coordinates (default: "z_location").
#' @param gene Character; column name for gene/feature identifiers (default: "feature_name").
#' @param qv Character; column name for quality values (default: "qv").
#' @param qv_threshold Numeric; minimum quality value to retain transcripts (default: 20).
#' @param n_neighbors Integer; number of nearest neighbors to connect (default: 10).
#' @param dist_threshold Numeric; maximum spatial distance for edge inclusion (default: 2).
#' @param same_label_ratio Numeric; log-probability ratio for same vs. different labels
#'   in the pairwise potential (default: 5).
#' @param ... Additional arguments passed to `potts_lbp_parallel_cpp()`.
#'
#' @return A list containing:
#'   \item{marginals}{Matrix of posterior marginal probabilities (transcripts x cell types).}
#'   \item{transcripts_df}{The input data frame with added `label` column containing
#'     assigned cell type labels.}
#'
#' @examples
#' \dontrun{
#' # Assuming transcripts_df and cell_signatures are prepared
#' result <- run_crf(
#'   transcripts_df,
#'   cell_signatures,
#'   qv_threshold = 20,
#'   n_neighbors = 10
#' )
#' head(result$transcripts_df)
#' }
#'
#' @export
run_crf = function(
    transcripts_df,
    cell_signatures,
    x = "x_location",
    y = "y_location",
    z = "z_location",
    gene = "feature_name",
    qv = "qv",
    qv_threshold = 20,
    n_neighbors = 10,
    dist_threshold = 2,
    same_label_ratio = 5,
    ...
) {

    # Filter transcripts by quality threshold
    transcripts_df = transcripts_df %>% filter(!!sym(qv) >= qv_threshold)

    # Compute Hilbert indices for spatial coherence ordering
    coords = transcripts_df[, c(x, y)] %>% as.matrix()
    index = hilbert_index(coords)
    transcripts_df = transcripts_df %>% mutate(hilbert = index) %>%
        arrange(.data$hilbert) %>%
        mutate(index = row_number())

    # Find k-nearest neighbors in 3D space
    transcript_adj = RANN::nn2(
        transcripts_df[, c(x, y, z)],
        k = n_neighbors + 1
    )

    # Extract edge list from neighbor indices and distances
    from = rep(seq_len(nrow(transcript_adj$nn.idx)), n_neighbors)
    to = as.numeric(transcript_adj$nn.idx[, -1])
    idx = which(as.numeric(transcript_adj$nn.dists[, -1]) <= dist_threshold)
    from = from[idx]
    to = to[idx]

    # Build sparse adjacency matrix with edge weights representing label agreement
    adj = Matrix::sparseMatrix(
        i = from,
        j = to,
        x = 1,
        dims = c(nrow(transcripts_df), nrow(transcripts_df))
    )
    adj@x = rep(log(same_label_ratio), length(adj@x))

    # Convert adjacency matrix to CSR-like graph format for LBP
    graph = build_potts_lbp_graph(adj)

    # Extract node potentials (log-likelihoods) from cell signatures
    node_potentials = log(cell_signatures[transcripts_df[[gene]], , drop = FALSE])

    # Run parallel loopy belief propagation for inference
    res = potts_lbp_parallel_cpp(
        graph$adj_ptr,
        graph$adj_idx,
        graph$rev_idx,
        graph$edge_weights,
        node_potentials,
        ...
    )

    # Extract maximum a posteriori (MAP) labels from marginals
    labels = apply(res$marginals, 1, which.max)
    transcripts_df = transcripts_df %>% mutate(label = factor(colnames(cell_signatures)[labels]))
    res$transcripts_df = transcripts_df
    colnames(res$marginals) = colnames(cell_signatures)

    return(res)
}
