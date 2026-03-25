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
    transcripts_df = transcripts_df %>% filter(!!sym(qv) >= qv_threshold)

    coords = transcripts_df[, c(x, y)] %>% as.matrix()
    index = hilbert_index(coords)
    transcripts_df = transcripts_df %>% mutate(hilbert = index) %>%
        arrange(hilbert) %>%
        mutate(index = row_number())

    transcript_adj = RANN::nn2(
        transcripts_df[, c(x, y, z)],
        k = n_neighbors + 1
    )
    from = rep(seq_len(nrow(transcript_adj$nn.idx)), n_neighbors)
    to = as.numeric(transcript_adj$nn.idx[, -1])
    idx = which(as.numeric(transcript_adj$nn.dists[, -1]) <= dist_threshold)
    from = from[idx]
    to = to[idx]

    adj = sparseMatrix::sparseMatrix(
        i = from,
        j = to,
        x = 1,
        dims = c(nrow(transcripts_df), nrow(transcripts_df))
    )
    adj@x = rep(log(same_label_ratio), length(adj@x))

    graph = build_potts_lbp_graph(adj)
    node_potentials = log(cell_signatures[transcripts_df[[gene]], , drop = FALSE])

    res = potts_lbp_parallel_cpp(
        graph$adj_ptr,
        graph$adj_idx,
        graph$rev_idx,
        graph$edge_weights,
        node_potentials,
        ...
    )

    labels = apply(res$marginals, 1, which.max)
    transcripts_df = transcripts_df %>% mutate(label = factor(colnames(cell_signatures)[labels]))
    res$transcripts_df = transcripts_df
    colnames(res$marginals) = colnames(cell_signatures)

    return(res)
}
