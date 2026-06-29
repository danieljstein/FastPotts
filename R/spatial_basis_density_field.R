build_active_edge_design <- function(basis_id) {
    n_active = ncol(basis_id)
    pairs = utils::combn(seq_len(n_active), 2L)

    pair_key = matrix(NA_character_, nrow = nrow(basis_id), ncol = ncol(pairs))
    edge_id = matrix(NA_integer_, nrow = nrow(basis_id), ncol = ncol(pairs))
    for (p in seq_len(ncol(pairs))) {
        a = pairs[1L, p]
        b = pairs[2L, p]
        pair_key[, p] = paste(
            pmin(basis_id[, a], basis_id[, b]),
            pmax(basis_id[, a], basis_id[, b]),
            sep = ":"
        )
    }

    unique_key = unique(as.vector(pair_key))
    edge_lookup = seq_along(unique_key)
    names(edge_lookup) = unique_key
    edge_id[] = unname(edge_lookup[pair_key])

    edge_pairs = do.call(rbind, strsplit(unique_key, ":", fixed = TRUE))
    storage.mode(edge_pairs) = "integer"
    colnames(edge_pairs) = c("from", "to")

    list(
        edge_id = edge_id - 1L,
        edge_pairs = edge_pairs,
        pair_from = as.integer(pairs[1L, ] - 1L),
        pair_to = as.integer(pairs[2L, ] - 1L)
    )
}

make_density_quadrature_grid <- function(coords, n, expansion = 0.02) {
    d = ncol(coords)
    mins = apply(coords, 2L, min)
    maxs = apply(coords, 2L, max)
    widths = pmax(maxs - mins, .Machine$double.eps)
    mins = mins - expansion * widths
    maxs = maxs + expansion * widths
    widths = maxs - mins

    n_axis = ceiling(n^(1 / d))
    axes = lapply(seq_len(d), function(j) seq(mins[j], maxs[j], length.out = n_axis))
    grid = as.matrix(expand.grid(axes))
    storage.mode(grid) = "double"
    quad_weight = rep(prod(widths) / nrow(grid), nrow(grid))

    list(coords = grid, weight = quad_weight, bounds = rbind(min = mins, max = maxs))
}

normalize_density_floor <- function(density_floor, cell_types) {
    if (length(density_floor) == 1L) {
        density_floor = rep(density_floor, length(cell_types))
    }
    density_floor = as.numeric(density_floor)
    if (length(density_floor) != length(cell_types)) {
        stop("density_floor must have length 1 or one value per posterior column.", call. = FALSE)
    }
    if (any(!is.finite(density_floor)) || any(density_floor < 0) || any(density_floor >= 1)) {
        stop("density_floor values must be in [0, 1).", call. = FALSE)
    }
    names(density_floor) = cell_types
    density_floor
}

#' Fit a spatial transcript density field with boundary attenuation
#'
#' Experimental prototype for fitting a posterior-weighted transcript density
#' field on the same triangular or BCC barycentric basis used by
#' [spatial_basis_segmentation()]. For each cell type, the modeled density is
#'
#' \deqn{
#' \rho_k(x) = \exp(\eta_k(x))
#' \left[f_k + (1 - f_k)\operatorname{sigmoid}(-g_k(x))\right],
#' }
#'
#' where `eta` is a vertex-linear baseline log-density field, and `g` combines
#' a vertex-linear field with shared edge-quadratic terms,
#'
#' \deqn{
#' g_k(x) = \sum_a \lambda_a(x) s_{a,k}
#' + \sum_{a<b} \lambda_a(x)\lambda_b(x) s_{ab,k}.
#' }
#'
#' The objective is a posterior-weighted inhomogeneous Poisson point-process
#' likelihood with a regular-grid quadrature approximation over the coordinate
#' bounding box.
#'
#' @param transcripts_df Transcript-level data frame.
#' @param posterior Numeric transcript-by-cell-type posterior probability
#'   matrix.
#' @param basis Character; `"2d"`/`"tri"` or `"3d"`/`"bcc"`.
#' @param s Positive spatial basis mesh size.
#' @param x,y,z Character coordinate column names.
#' @param origin Optional lattice origin.
#' @param quadrature_n Approximate number of quadrature points over the
#'   coordinate bounding box.
#' @param quadrature_expansion Fractional bounding-box expansion used for
#'   quadrature.
#' @param density_floor Scalar or vector floor multiplier in `[0, 1)`. A value
#'   of `0.05` means boundary attenuation can reduce density to 5 percent of
#'   the baseline.
#' @param lambda_eta Smoothness penalty for the baseline log-density field.
#' @param lambda_s_smooth Smoothness penalty for the vertex-linear boundary
#'   logit field.
#' @param lambda_s_prior L2 shrinkage strength for vertex boundary logits
#'   toward `s_prior_mean`.
#' @param s_prior_mean Prior mean for vertex boundary logits. Negative values
#'   favor no boundary attenuation.
#' @param lambda_edge L2 shrinkage strength for edge-quadratic boundary terms.
#' @param maxit Maximum L-BFGS iterations.
#' @param reltol Relative convergence tolerance.
#' @param n_threads Integer number of OpenMP threads. If `NULL`, uses runtime
#'   default.
#' @param show_progress Logical; print progress messages.
#'
#' @return A list with fitted parameters, transcript-level predictions,
#'   quadrature design, basis design, optimizer result, and effective
#'   parameters.
#' @export
spatial_basis_density_field <- function(
    transcripts_df,
    posterior,
    basis = c("3d", "2d", "tri", "bcc"),
    s,
    x = "x_location",
    y = "y_location",
    z = "z_location",
    origin = NULL,
    quadrature_n = 10000L,
    quadrature_expansion = 0.02,
    density_floor = 0.05,
    lambda_eta = 0.1,
    lambda_s_smooth = 0.1,
    lambda_s_prior = 0.1,
    s_prior_mean = -4,
    lambda_edge = 0.25,
    maxit = 100L,
    reltol = 1e-6,
    n_threads = NULL,
    show_progress = TRUE
) {
    basis = normalize_spatial_basis(basis)
    lattice_basis = if (basis == "2d") "tri" else "bcc"
    d = if (basis == "2d") 2L else 3L
    coord_cols = if (basis == "2d") c(x, y) else c(x, y, z)
    missing_cols = setdiff(coord_cols, colnames(transcripts_df))
    if (length(missing_cols) > 0L) {
        stop("Missing coordinate column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    if (is.null(origin)) {
        origin = rep(0, d)
    }
    if (length(origin) != d || any(!is.finite(origin))) {
        stop("origin must be a finite numeric vector with length matching the selected basis.", call. = FALSE)
    }
    if (length(s) != 1L || !is.finite(s) || s <= 0) {
        stop("s must be a positive finite scalar.", call. = FALSE)
    }
    if (quadrature_n < 1L) {
        stop("quadrature_n must be positive.", call. = FALSE)
    }
    if (!is.finite(quadrature_expansion) || quadrature_expansion < 0) {
        stop("quadrature_expansion must be non-negative and finite.", call. = FALSE)
    }
    for (arg in c("lambda_eta", "lambda_s_smooth", "lambda_s_prior", "lambda_edge")) {
        value = get(arg)
        if (length(value) != 1L || !is.finite(value) || value < 0) {
            stop(arg, " must be a non-negative finite scalar.", call. = FALSE)
        }
    }
    if (length(s_prior_mean) != 1L || !is.finite(s_prior_mean)) {
        stop("s_prior_mean must be a finite scalar.", call. = FALSE)
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

    posterior = normalize_posterior_matrix(posterior, nrow(transcripts_df))
    cell_types = colnames(posterior)
    if (is.null(cell_types)) {
        cell_types = paste0("type", seq_len(ncol(posterior)))
        colnames(posterior) = cell_types
    }
    density_floor = normalize_density_floor(density_floor, cell_types)

    coords = as.matrix(transcripts_df[, coord_cols, drop = FALSE])
    storage.mode(coords) = "double"
    if (any(!is.finite(coords))) {
        stop("Coordinate columns must contain finite values.", call. = FALSE)
    }

    if (show_progress) {
        message("Building quadrature grid...")
    }
    quadrature = make_density_quadrature_grid(coords, quadrature_n, quadrature_expansion)
    all_coords = rbind(coords, quadrature$coords)

    if (show_progress) {
        message("Computing spatial basis interpolation...")
    }
    bary = if (basis == "2d") {
        tri_barycentric(all_coords, s = s, origin = origin, n_threads = basis_n_threads)
    } else {
        bcc_barycentric(all_coords, s = s, origin = origin, n_threads = basis_n_threads)
    }
    design = basis_design_from_barycentric(bary)
    n_obs = nrow(coords)
    obs_idx = seq_len(n_obs)
    quad_idx = seq.int(n_obs + 1L, nrow(all_coords))

    if (show_progress) {
        message("Building basis edge design...")
    }
    basis_edges = build_lattice_neighbor_edges(design$basis_lattice, basis = lattice_basis, s = s)
    active_edge = build_active_edge_design(design$basis_id)

    obs_basis_id = design$basis_id[obs_idx, , drop = FALSE] - 1L
    quad_basis_id = design$basis_id[quad_idx, , drop = FALSE] - 1L
    obs_basis_weight = design$basis_weight[obs_idx, , drop = FALSE]
    quad_basis_weight = design$basis_weight[quad_idx, , drop = FALSE]
    obs_edge_id = active_edge$edge_id[obs_idx, , drop = FALSE]
    quad_edge_id = active_edge$edge_id[quad_idx, , drop = FALSE]

    M = nrow(design$basis_lattice)
    E = nrow(active_edge$edge_pairs)
    K = ncol(posterior)
    par0 = numeric((2L * M + E) * K)
    eta_offset = 0L
    s_offset = M * K
    edge_offset = 2L * M * K

    # Initialize the broad density near the observed posterior mass per volume.
    volume = sum(quadrature$weight)
    type_mass = colSums(posterior)
    eta0 = log(pmax(type_mass / volume, 1e-8))
    for (k in seq_len(K)) {
        par0[eta_offset + seq_len(M) + M * (k - 1L)] = eta0[k]
        par0[s_offset + seq_len(M) + M * (k - 1L)] = s_prior_mean
    }

    objective = function(par) {
        spatial_density_objective_cpp(
            par = par,
            obs_basis_id = obs_basis_id,
            obs_basis_weight = obs_basis_weight,
            obs_edge_id = obs_edge_id,
            obs_weight = posterior,
            quad_basis_id = quad_basis_id,
            quad_basis_weight = quad_basis_weight,
            quad_edge_id = quad_edge_id,
            quad_weight = quadrature$weight,
            pair_from = active_edge$pair_from,
            pair_to = active_edge$pair_to,
            smooth_edge_from = as.integer(basis_edges[, "from"] - 1L),
            smooth_edge_to = as.integer(basis_edges[, "to"] - 1L),
            smooth_edge_distance = as.numeric(basis_edges[, "distance"]),
            density_floor = as.numeric(density_floor),
            lambda_eta = lambda_eta,
            lambda_s_smooth = lambda_s_smooth,
            lambda_s_prior = lambda_s_prior,
            s_prior_mean = s_prior_mean,
            lambda_edge = lambda_edge,
            n_basis = M,
            n_edge_coef = E,
            n_cell_types = K,
            n_threads = n_threads
        )
    }

    if (show_progress) {
        message("Optimizing spatial density field...")
    }
    opt = stats::optim(
        par = par0,
        fn = function(par) objective(par)$value,
        gr = function(par) objective(par)$gradient,
        method = "L-BFGS-B",
        control = list(maxit = as.integer(maxit), factr = reltol / .Machine$double.eps)
    )
    warn_spatial_basis_optim_status(opt, maxit)

    pred = spatial_density_predict_cpp(
        par = opt$par,
        basis_id = obs_basis_id,
        basis_weight = obs_basis_weight,
        edge_id = obs_edge_id,
        pair_from = active_edge$pair_from,
        pair_to = active_edge$pair_to,
        density_floor = as.numeric(density_floor),
        n_basis = M,
        n_edge_coef = E,
        n_cell_types = K,
        n_threads = n_threads
    )
    for (name in names(pred)) {
        colnames(pred[[name]]) = cell_types
    }

    eta = matrix(opt$par[eta_offset + seq_len(M * K)], nrow = M, ncol = K)
    s_node = matrix(opt$par[s_offset + seq_len(M * K)], nrow = M, ncol = K)
    s_edge = matrix(opt$par[edge_offset + seq_len(E * K)], nrow = E, ncol = K)
    colnames(eta) = colnames(s_node) = colnames(s_edge) = cell_types

    list(
        density = pred$density,
        eta = pred$eta,
        boundary_logit = pred$g,
        attenuation = pred$attenuation,
        eta_basis = eta,
        boundary_node_basis = s_node,
        boundary_edge_basis = s_edge,
        density_floor = density_floor,
        basis_points = design$basis_points,
        basis_lattice = design$basis_lattice,
        basis_edges = basis_edges,
        quadrature = quadrature,
        optim = opt,
        parameters = list(
            basis = basis,
            s = s,
            origin = origin,
            quadrature_n = quadrature_n,
            quadrature_expansion = quadrature_expansion,
            lambda_eta = lambda_eta,
            lambda_s_smooth = lambda_s_smooth,
            lambda_s_prior = lambda_s_prior,
            s_prior_mean = s_prior_mean,
            lambda_edge = lambda_edge,
            maxit = maxit,
            reltol = reltol
        )
    )
}
