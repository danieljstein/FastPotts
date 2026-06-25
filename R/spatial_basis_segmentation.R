basis_design_from_barycentric <- function(basis) {
    lattice = basis$lattice
    dims = dim(lattice)
    n = dims[1]
    n_active = dims[2]
    d = dims[3]

    flat_lattice = do.call(
        cbind,
        lapply(seq_len(d), function(j) as.vector(lattice[, , j]))
    )
    key = do.call(paste, c(as.data.frame(flat_lattice), sep = ":"))
    unique_key = unique(key)
    basis_id = matrix(match(key, unique_key), nrow = n, ncol = n_active)

    basis_lattice = flat_lattice[match(unique_key, key), , drop = FALSE]
    basis_points = do.call(
        cbind,
        lapply(seq_len(d), function(j) as.vector(basis$points[, , j])[match(unique_key, key)])
    )

    list(
        basis_id = basis_id,
        basis_weight = basis$weights,
        basis_lattice = basis_lattice,
        basis_points = basis_points
    )
}

build_lattice_neighbor_edges <- function(basis_lattice, basis, s) {
    basis = match.arg(basis, c("tri", "bcc"))
    key = do.call(paste, c(as.data.frame(basis_lattice), sep = ":"))
    key_to_id = seq_along(key)
    names(key_to_id) = key

    if (basis == "tri") {
        offsets = rbind(
            c(1L, 0L),
            c(0L, 1L),
            c(1L, -1L)
        )
        offset_lengths = rep(as.numeric(s), nrow(offsets))
    } else {
        offsets = rbind(
            c(2L, 0L, 0L),
            c(0L, 2L, 0L),
            c(0L, 0L, 2L),
            as.matrix(expand.grid(c(-1L, 1L), c(-1L, 1L), c(-1L, 1L)))
        )
        offset_lengths = sqrt(rowSums((offsets * as.numeric(s) / 2)^2))
    }

    from = integer()
    to = integer()
    distance = numeric()

    for (a in seq_len(nrow(offsets))) {
        neighbor = sweep(basis_lattice, 2, offsets[a, ], "+")
        neighbor_key = do.call(paste, c(as.data.frame(neighbor), sep = ":"))
        neighbor_id = unname(key_to_id[neighbor_key])
        keep = !is.na(neighbor_id)

        from = c(from, which(keep))
        to = c(to, neighbor_id[keep])
        distance = c(distance, rep(offset_lengths[a], sum(keep)))
    }

    keep = from < to
    cbind(from = from[keep], to = to[keep], distance = distance[keep])
}

normalize_spatial_basis <- function(basis) {
    basis = match.arg(tolower(basis), c("2d", "3d", "tri", "bcc"))
    if (basis == "tri") {
        return("2d")
    }
    if (basis == "bcc") {
        return("3d")
    }
    basis
}

update_signatures_from_posteriors <- function(
    gene_index,
    posterior,
    reference_signatures,
    current_signatures,
    prior_strength,
    update_rate,
    min_posterior,
    signature_floor
) {
    G = nrow(reference_signatures)
    K = ncol(reference_signatures)
    weights = posterior

    if (min_posterior > 0) {
        weights[weights < min_posterior] = 0
    }

    counts = matrix(0, nrow = G, ncol = K)
    for (k in seq_len(K)) {
        grouped = rowsum(weights[, k], group = gene_index, reorder = FALSE)
        counts[as.integer(rownames(grouped)), k] = grouped[, 1]
    }

    posterior_mean = sweep(
        counts + prior_strength * reference_signatures,
        2,
        colSums(counts) + prior_strength,
        "/"
    )
    updated = (1 - update_rate) * current_signatures + update_rate * posterior_mean
    updated = pmax(updated, signature_floor)
    updated = sweep(updated, 2, colSums(updated), "/")
    dimnames(updated) = dimnames(reference_signatures)

    list(
        signatures = updated,
        effective_counts = colSums(weights),
        max_abs_change = apply(abs(updated - current_signatures), 2, max),
        posterior_mean = posterior_mean
    )
}

warn_spatial_basis_optim_status <- function(opt, maxit) {
    if (is.null(opt$convergence) || opt$convergence == 0L) {
        return(invisible(NULL))
    }

    if (opt$convergence == 1L) {
        warning(
            "Optimization reached maxit = ",
            as.integer(maxit),
            " before satisfying the L-BFGS-B convergence criterion. ",
            "The returned fit is the best iterate found; consider increasing maxit ",
            "or checking fit$optim$convergence and fit$optim$value.",
            call. = FALSE
        )
        return(invisible(NULL))
    }

    if (opt$convergence == 51L) {
        warning(
            "L-BFGS-B reported a warning",
            if (!is.null(opt$message) && nzchar(opt$message)) paste0(": ", opt$message) else ".",
            call. = FALSE
        )
        return(invisible(NULL))
    }

    if (opt$convergence == 52L) {
        warning(
            "L-BFGS-B reported an error",
            if (!is.null(opt$message) && nzchar(opt$message)) paste0(": ", opt$message) else ".",
            call. = FALSE
        )
        return(invisible(NULL))
    }

    warning(
        "Optimization did not converge; optim() returned convergence code ",
        opt$convergence,
        if (!is.null(opt$message) && nzchar(opt$message)) paste0(" with message: ", opt$message) else ".",
        call. = FALSE
    )
    invisible(NULL)
}

#' Segment transcripts with a continuous spatial basis model
#'
#' Fits a continuous spatial cell-type field on either a 2D triangular lattice
#' or a 3D body-centered cubic (BCC) lattice. Each transcript has a spatial
#' prior over cell types,
#' `softmax(f_1(x_i), ..., f_K(x_i))`, where each `f_k(x_i)` is interpolated
#' from nearby lattice basis coefficients using barycentric weights. The
#' observed transcript gene then contributes the cell-type signature likelihood.
#'
#' This implementation performs MAP estimation of the spatial field with
#' optional conservative refinement of the cell-type signatures. The final
#' transcript posterior is proportional to the fitted spatial prior multiplied
#' by the corresponding gene signature probability.
#'
#' @param transcripts_df A data frame containing transcript-level data.
#' @param cell_signatures A numeric matrix with genes in rows and cell types in
#'   columns. Values are treated as gene emission probabilities for each cell
#'   type.
#' @param basis Character; either `"2d"` for a 2D triangular lattice or `"3d"`
#'   for a 3D BCC lattice. The older values `"tri"` and `"bcc"` are accepted as
#'   aliases.
#' @param s Positive numeric mesh size. For `basis = "2d"`, this is the side
#'   length of each equilateral Delaunay triangle. For `basis = "3d"`, this is
#'   the distance between same-parity BCC lattice points along a coordinate
#'   axis, matching the longest Delaunay tetrahedron edge length.
#' @param x Character; column name for x-coordinates.
#' @param y Character; column name for y-coordinates.
#' @param z Character; column name for z-coordinates. Used only when
#'   `basis = "3d"`.
#' @param gene Character; column name for gene identifiers.
#' @param qv Character; column name for quality values.
#' @param is_gene Character; column name indicating whether a transcript should
#'   be treated as a gene. If absent, no gene-only filtering is applied.
#' @param qv_threshold Numeric; minimum quality value to retain transcripts.
#' @param origin Numeric lattice origin. Defaults to zeros with length matching
#'   the selected basis dimension.
#' @param lambda Non-negative smoothing strength between neighboring basis
#'   coefficients. The transcript likelihood is averaged over transcripts and
#'   the spatial penalty is averaged over graph edges, so `lambda` is on an
#'   approximate per-transcript, per-edge spatial-slope scale.
#' @param regularization Character; one of `"quadratic"`, `"huber"`, or
#'   `"bounded"`. All options are applied to the spatial logit slope between
#'   neighboring basis points, `(w_m - w_n) / edge_distance`. Quadratic
#'   smoothing penalizes squared slope, Huber smoothing is quadratic near zero
#'   and linear past `delta`, and bounded smoothing has quadratic small-slope
#'   behavior with asymptotic bound `sigma^2`.
#' @param delta Positive Huber transition scale in logit units per spatial
#'   unit. Used only when `regularization = "huber"`.
#' @param sigma Positive bounded-penalty slope scale in logit units per spatial
#'   unit. Used only when `regularization = "bounded"`; the per-component
#'   asymptotic penalty is `sigma^2`.
#' @param purity Character; one of `"none"`, `"entropy"`, or `"gini"`.
#'   Entropy and Gini purity penalties encourage each basis point's cell-type
#'   prior to concentrate on fewer cell types.
#' @param purity_lambda Non-negative strength of the basis-point purity penalty.
#' @param signature_floor Positive floor applied to signatures before
#'   log-transforming.
#' @param normalize_signatures Logical; if `TRUE`, normalize each cell-type
#'   signature column after applying `signature_floor`.
#' @param refine_signatures Logical; if `TRUE`, alternates spatial field
#'   fitting with conservative cell-type signature updates from posterior
#'   transcript assignments.
#' @param signature_update_iters Integer number of signature update steps. Each
#'   update is followed by a refit of the spatial field.
#' @param signature_prior_strength Non-negative Dirichlet pseudo-count strength
#'   for anchoring each refined signature to the input reference signature. If
#'   `NULL`, defaults to the number of genes in `cell_signatures`.
#' @param signature_update_rate Numeric in `[0, 1]`; damping rate for each
#'   signature update.
#' @param signature_min_posterior Numeric in `[0, 1]`; if positive, only
#'   posterior assignment weights at least this large contribute to signature
#'   soft counts.
#' @param maxit Integer maximum number of L-BFGS iterations.
#' @param reltol Approximate relative convergence tolerance. For the
#'   `"L-BFGS-B"` optimizer this is converted to `factr = reltol /
#'   .Machine$double.eps`.
#' @param n_threads Integer number of OpenMP threads for objective/gradient and
#'   prediction calculations. If `NULL`, uses the OpenMP runtime default.
#' @param show_progress Logical; if `TRUE`, prints progress messages.
#'
#' @return A list with components:
#' \describe{
#'   \item{marginals}{Matrix of posterior transcript probabilities.}
#'   \item{spatial_prior}{Matrix of fitted spatial prior probabilities.}
#'   \item{logits}{Matrix of fitted spatial logits at transcript locations.}
#'   \item{basis_weights}{Matrix of fitted basis coefficients, one row per
#'     lattice basis point and one column per cell type. The final cell type is
#'     the reference class with coefficient zero.}
#'   \item{basis_points}{Matrix of lattice basis point coordinates.}
#'   \item{basis_lattice}{Matrix of integer lattice coordinates.}
#'   \item{basis_edges}{Matrix of neighboring basis point indices and physical
#'     edge distances.}
#'   \item{transcripts_df}{Filtered input data with added `label` column.}
#'   \item{optim}{The [stats::optim()] result.}
#'   \item{optim_history}{List of optimizer results, one per spatial field fit.}
#'   \item{cell_signatures_initial}{Input signatures after filtering,
#'     flooring, and optional normalization.}
#'   \item{cell_signatures}{Final signatures used for the returned posterior.}
#'   \item{signature_history}{List of signatures after each refinement step.}
#'   \item{signature_update_history}{List of per-update diagnostics.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- spatial_basis_segmentation(
#'   transcripts_df,
#'   cell_signatures,
#'   basis = "2d",
#'   s = 2
#' )
#' fit$marginals
#' }
#'
#' @export
spatial_basis_segmentation <- function(
    transcripts_df,
    cell_signatures,
    basis = c("2d", "3d", "tri", "bcc"),
    s,
    x = "x_location",
    y = "y_location",
    z = "z_location",
    gene = "feature_name",
    qv = "qv",
    is_gene = "is_gene",
    qv_threshold = 20,
    origin = NULL,
    lambda = 1,
    regularization = c("quadratic", "huber", "bounded"),
    delta = 1,
    sigma = 1,
    purity = c("none", "entropy", "gini"),
    purity_lambda = 0,
    signature_floor = 1e-12,
    normalize_signatures = TRUE,
    refine_signatures = FALSE,
    signature_update_iters = 1L,
    signature_prior_strength = NULL,
    signature_update_rate = 0.25,
    signature_min_posterior = 0,
    maxit = 100L,
    reltol = 1e-6,
    n_threads = NULL,
    show_progress = TRUE
) {
    basis = normalize_spatial_basis(basis)
    lattice_basis = if (basis == "2d") "tri" else "bcc"
    regularization = match.arg(regularization)
    purity = match.arg(purity)
    d = if (basis == "2d") 2L else 3L

    if (is.null(origin)) {
        origin = rep(0, d)
    }
    if (length(origin) != d || any(!is.finite(origin))) {
        stop("origin must be a finite numeric vector with length matching the selected basis.")
    }
    if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
        stop("lambda must be a non-negative finite number.")
    }
    if (length(delta) != 1L || !is.finite(delta) || delta <= 0) {
        stop("delta must be a positive finite number.")
    }
    if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
        stop("sigma must be a positive finite number.")
    }
    if (length(purity_lambda) != 1L || !is.finite(purity_lambda) || purity_lambda < 0) {
        stop("purity_lambda must be a non-negative finite number.")
    }
    if (length(signature_floor) != 1L || !is.finite(signature_floor) || signature_floor <= 0) {
        stop("signature_floor must be a positive finite number.")
    }
    if (!is.logical(refine_signatures) || length(refine_signatures) != 1L || is.na(refine_signatures)) {
        stop("refine_signatures must be TRUE or FALSE.")
    }
    if (
        length(signature_update_iters) != 1L ||
        !is.finite(signature_update_iters) ||
        signature_update_iters < 0 ||
        signature_update_iters != as.integer(signature_update_iters)
    ) {
        stop("signature_update_iters must be a non-negative integer.")
    }
    signature_update_iters = as.integer(signature_update_iters)
    if (length(signature_update_rate) != 1L || !is.finite(signature_update_rate) || signature_update_rate < 0 || signature_update_rate > 1) {
        stop("signature_update_rate must be a finite number in [0, 1].")
    }
    if (length(signature_min_posterior) != 1L || !is.finite(signature_min_posterior) || signature_min_posterior < 0 || signature_min_posterior > 1) {
        stop("signature_min_posterior must be a finite number in [0, 1].")
    }
    if (is.null(n_threads)) {
        n_threads = 0L
    } else if (
        length(n_threads) != 1L ||
        !is.finite(n_threads) ||
        n_threads < 1 ||
        n_threads != as.integer(n_threads)
    ) {
        stop("n_threads must be NULL or a positive integer.")
    } else {
        n_threads = as.integer(n_threads)
    }
    basis_n_threads = if (n_threads == 0L) NULL else n_threads

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
    if (sum(keep_signature) < nrow(df)) {
        warning(
            "Only ",
            sum(keep_signature),
            " out of ",
            nrow(df),
            " transcripts have genes found in cell_signatures. Filtering to these transcripts."
        )
    }
    df = df[keep_signature, , drop = FALSE]

    if (nrow(df) == 0L) {
        stop("No transcripts remain after filtering.")
    }

    signatures = as.matrix(cell_signatures)
    storage.mode(signatures) = "double"
    if (any(!is.finite(signatures)) || any(signatures < 0)) {
        stop("cell_signatures must contain finite non-negative values.")
    }
    signatures = pmax(signatures, signature_floor)
    if (isTRUE(normalize_signatures)) {
        signatures = sweep(signatures, 2, colSums(signatures), "/")
    }
    reference_signatures = signatures
    if (is.null(signature_prior_strength)) {
        signature_prior_strength = nrow(signatures)
    } else if (length(signature_prior_strength) != 1L || !is.finite(signature_prior_strength) || signature_prior_strength < 0) {
        stop("signature_prior_strength must be NULL or a non-negative finite number.")
    }

    coord_cols = if (basis == "2d") c(x, y) else c(x, y, z)
    coords = as.matrix(df[, coord_cols, drop = FALSE])

    if (show_progress) {
        message("Computing spatial basis interpolation...")
    }
    bary = if (basis == "2d") {
        tri_barycentric(coords, s = s, origin = origin, n_threads = basis_n_threads)
    } else {
        bcc_barycentric(coords, s = s, origin = origin, n_threads = basis_n_threads)
    }
    design = basis_design_from_barycentric(bary)

    if (show_progress) {
        message("Building lattice neighbor graph...")
    }
    basis_edges = build_lattice_neighbor_edges(design$basis_lattice, basis = lattice_basis, s = s)

    gene_index = match(df[[gene]], rownames(signatures))
    K = ncol(signatures)
    M = nrow(design$basis_lattice)
    if (K < 2L) {
        stop("cell_signatures must contain at least two cell types.")
    }

    par0 = numeric(M * (K - 1L))
    basis_id0 = design$basis_id - 1L
    edge_from0 = as.integer(basis_edges[, "from"] - 1L)
    edge_to0 = as.integer(basis_edges[, "to"] - 1L)
    edge_distance = as.numeric(basis_edges[, "distance"])
    regularization_id = match(regularization, c("quadratic", "huber", "bounded")) - 1L
    purity_id = match(purity, c("none", "entropy", "gini")) - 1L

    fit_spatial_field = function(par_start, current_signatures, fit_iter) {
        log_signature = log(current_signatures)
        objective = function(par) {
            spatial_basis_objective_cpp(
                par = par,
                basis_id = basis_id0,
                basis_weight = design$basis_weight,
                gene_index = as.integer(gene_index - 1L),
                log_signature = log_signature,
                edge_from = edge_from0,
                edge_to = edge_to0,
                edge_distance = edge_distance,
                lambda = lambda,
                regularization = regularization_id,
                delta = delta,
                sigma = sigma,
                purity = purity_id,
                purity_lambda = purity_lambda,
                n_threads = n_threads,
                n_basis = M,
                n_cell_types = K
            )
        }

        if (show_progress) {
            message("Optimizing continuous spatial field", if (fit_iter > 1L) paste0(" (fit ", fit_iter, ")") else "", "...")
        }
        opt = stats::optim(
            par = par_start,
            fn = function(par) objective(par)$value,
            gr = function(par) objective(par)$gradient,
            method = "L-BFGS-B",
            control = list(
                maxit = as.integer(maxit),
                factr = reltol / .Machine$double.eps
            )
        )
        warn_spatial_basis_optim_status(opt, maxit)

        pred = spatial_basis_predict_cpp(
            par = opt$par,
            basis_id = basis_id0,
            basis_weight = design$basis_weight,
            gene_index = as.integer(gene_index - 1L),
            log_signature = log_signature,
            n_basis = M,
            n_threads = n_threads,
            n_cell_types = K
        )

        list(opt = opt, pred = pred)
    }

    n_updates = if (isTRUE(refine_signatures)) signature_update_iters else 0L
    n_fits = n_updates + 1L
    current_signatures = signatures
    par_start = par0
    signature_history = list(current_signatures)
    optim_history = vector("list", n_fits)
    signature_update_history = vector("list", n_updates)

    for (fit_iter in seq_len(n_fits)) {
        fit = fit_spatial_field(par_start, current_signatures, fit_iter)
        opt = fit$opt
        pred = fit$pred
        optim_history[[fit_iter]] = opt

        if (fit_iter <= n_updates) {
            if (show_progress) {
                message("Updating cell type signatures...")
            }
            update = update_signatures_from_posteriors(
                gene_index = gene_index,
                posterior = pred$posterior,
                reference_signatures = reference_signatures,
                current_signatures = current_signatures,
                prior_strength = signature_prior_strength,
                update_rate = signature_update_rate,
                min_posterior = signature_min_posterior,
                signature_floor = signature_floor
            )
            current_signatures = update$signatures
            signature_history[[fit_iter + 1L]] = current_signatures
            signature_update_history[[fit_iter]] = update[c("effective_counts", "max_abs_change")]
            par_start = opt$par
        }
    }

    colnames(pred$posterior) = colnames(current_signatures)
    colnames(pred$prior) = colnames(current_signatures)
    colnames(pred$logits) = colnames(current_signatures)

    basis_weights = matrix(0, nrow = M, ncol = K)
    basis_weights[, seq_len(K - 1L)] = matrix(opt$par, nrow = M, ncol = K - 1L)
    colnames(basis_weights) = colnames(current_signatures)

    labels = max.col(pred$posterior, ties.method = "first")
    df$label = factor(colnames(current_signatures)[labels], levels = colnames(current_signatures))

    list(
        marginals = pred$posterior,
        spatial_prior = pred$prior,
        logits = pred$logits,
        basis_weights = basis_weights,
        basis_points = design$basis_points,
        basis_lattice = design$basis_lattice,
        basis_edges = basis_edges,
        transcripts_df = df,
        optim = opt,
        optim_history = optim_history,
        cell_signatures_initial = reference_signatures,
        cell_signatures = current_signatures,
        signature_history = signature_history,
        signature_update_history = signature_update_history
    )
}
