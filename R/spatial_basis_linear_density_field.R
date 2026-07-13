#' Fit a piecewise-linear positive spatial transcript density field
#'
#' Fits an inhomogeneous Poisson point-process density model on the same
#' triangular or BCC barycentric basis used by [spatial_basis_segmentation()].
#' Unlike [spatial_basis_log_density_field()], this model interpolates the
#' density/intensity itself:
#'
#' \deqn{
#' \rho_a = \exp(\alpha_a), \quad
#' \rho(x) = \sum_a \lambda_a(x) \rho_a.
#' }
#'
#' The integral term is exact over occupied simplexes because the integral of a
#' piecewise-linear basis function is the simplex volume divided by the number
#' of simplex vertices. Smoothness penalties are applied to the log-density
#' coefficients `alpha`, not to the linear intensities.
#'
#' @inheritParams spatial_basis_log_density_field
#'
#' @return A list with optional transcript-level density outputs, fitted
#'   basis-level density coefficients, exact simplex domain, basis design,
#'   transcript interpolation design, optimizer diagnostics, and effective
#'   parameters.
#' @export
spatial_basis_linear_density_field <- function(
    transcripts_df,
    posterior = NULL,
    return_density = FALSE,
    return_total_density = FALSE,
    return_eta = FALSE,
    store_posterior = FALSE,
    basis = c("3d", "2d", "tri", "bcc"),
    s,
    x = "x_location",
    y = "y_location",
    z = "z_location",
    origin = NULL,
    lambda = 0.1,
    regularization = c("bounded", "quadratic", "huber"),
    delta = 1,
    sigma = 1,
    lambda_laplacian = 0,
    domain_expansion_steps = 0L,
    domain_expansion_axes = c("xyz", "xy"),
    maxit = 100L,
    reltol = 1e-6,
    n_threads = NULL,
    show_progress = TRUE,
    optim = c("summary", "full")
) {
    basis = normalize_spatial_basis(basis)
    regularization = match.arg(regularization)
    domain_expansion_axes = match.arg(domain_expansion_axes)
    optim = match.arg(optim)
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
    if (!is.logical(return_density) || length(return_density) != 1L || is.na(return_density)) {
        stop("return_density must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(return_total_density) || length(return_total_density) != 1L || is.na(return_total_density)) {
        stop("return_total_density must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(return_eta) || length(return_eta) != 1L || is.na(return_eta)) {
        stop("return_eta must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(store_posterior) || length(store_posterior) != 1L || is.na(store_posterior)) {
        stop("store_posterior must be TRUE or FALSE.", call. = FALSE)
    }
    if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
        stop("lambda must be a non-negative finite scalar.", call. = FALSE)
    }
    if (length(delta) != 1L || !is.finite(delta) || delta <= 0) {
        stop("delta must be a positive finite scalar.", call. = FALSE)
    }
    if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
        stop("sigma must be a positive finite scalar.", call. = FALSE)
    }
    if (length(lambda_laplacian) != 1L || !is.finite(lambda_laplacian) || lambda_laplacian < 0) {
        stop("lambda_laplacian must be a non-negative finite scalar.", call. = FALSE)
    }
    if (
        length(domain_expansion_steps) != 1L ||
        !is.finite(domain_expansion_steps) ||
        domain_expansion_steps < 0 ||
        domain_expansion_steps != as.integer(domain_expansion_steps)
    ) {
        stop("domain_expansion_steps must be a non-negative integer.", call. = FALSE)
    }
    domain_expansion_steps = as.integer(domain_expansion_steps)
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

    coords = as.matrix(transcripts_df[, coord_cols, drop = FALSE])
    storage.mode(coords) = "double"
    if (any(!is.finite(coords))) {
        stop("Coordinate columns must contain finite values.", call. = FALSE)
    }

    if (!is.null(posterior)) {
        posterior = normalize_posterior_matrix(posterior, nrow(transcripts_df))
        cell_types = colnames(posterior)
        if (is.null(cell_types)) {
            cell_types = paste0("type", seq_len(ncol(posterior)))
            colnames(posterior) = cell_types
        }
    } else {
        cell_types = character()
    }
    if (return_density && is.null(posterior)) {
        stop("posterior must be supplied when return_density = TRUE.", call. = FALSE)
    }

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
        message("Building simplex integration domain...")
    }
    domain = make_log_density_domain_simplex(
        design = design,
        basis = lattice_basis,
        s = s,
        origin = origin,
        domain_expansion_steps = domain_expansion_steps,
        domain_expansion_axes = domain_expansion_axes
    )
    basis_lattice = domain$basis_lattice
    basis_points = domain$basis_points

    if (show_progress) {
        message("Building basis edge design...")
    }
    basis_edges = build_lattice_neighbor_edges(basis_lattice, basis = lattice_basis, s = s)
    lambda_laplacian_effective = lambda_laplacian / (s^4)
    edge_from = as.integer(basis_edges[, "from"] - 1L)
    edge_to = as.integer(basis_edges[, "to"] - 1L)
    edge_distance = as.numeric(basis_edges[, "distance"])

    obs_basis_id = design$basis_id - 1L
    obs_basis_weight = design$basis_weight
    simplex_basis_id = domain$basis_id - 1L

    volume = if (length(domain$volume) == 1L) {
        domain$volume * domain$n_simplex
    } else {
        sum(domain$volume)
    }
    par0 = rep(log(pmax(nrow(coords) / volume, 1e-8)), nrow(basis_lattice))
    regularization_id = match(regularization, c("quadratic", "huber", "bounded")) - 1L

    objective = function(par) {
        spatial_linear_density_objective_simplex_cpp(
            par = par,
            obs_basis_id = obs_basis_id,
            obs_basis_weight = obs_basis_weight,
            simplex_basis_id = simplex_basis_id,
            simplex_volume = domain$volume,
            edge_from = edge_from,
            edge_to = edge_to,
            edge_distance = edge_distance,
            lambda = lambda,
            regularization = regularization_id,
            delta = delta,
            sigma = sigma,
            lambda_laplacian = lambda_laplacian_effective,
            n_threads = n_threads
        )
    }

    if (show_progress) {
        message("Optimizing linear-density field...")
    }
    cached = make_cached_spatial_basis_objective(objective)
    opt = stats::optim(
        par = par0,
        fn = cached$fn,
        gr = cached$gr,
        method = "L-BFGS-B",
        control = list(maxit = as.integer(maxit), factr = reltol / .Machine$double.eps)
    )
    attr(opt, "objective_evaluations") = cached$n_eval()
    attr(opt, "objective_cache_hits") = cached$n_hit()
    warn_spatial_basis_optim_status(opt, maxit)

    need_transcript_prediction = return_density || return_total_density || return_eta
    pred = if (need_transcript_prediction) {
        spatial_linear_density_predict_cpp(
            par = opt$par,
            basis_id = obs_basis_id,
            basis_weight = obs_basis_weight,
            n_threads = n_threads
        )
    } else {
        NULL
    }

    out = list(
        density = NULL,
        total_density = NULL,
        eta = NULL,
        posterior = NULL,
        total_density_basis = as.numeric(exp(opt$par)),
        eta_basis = opt$par,
        basis_points = basis_points,
        basis_lattice = basis_lattice,
        basis_edges = basis_edges,
        transcript_basis_id = design$basis_id,
        transcript_basis_weight = design$basis_weight,
        domain = domain,
        quadrature = NULL,
        optim = if (optim == "full") opt else summarize_spatial_basis_optim(opt),
        parameters = list(
            mode = "linear_density",
            basis = basis,
            s = s,
            origin = origin,
            integration = "exact_linear_simplex",
            return_density = return_density,
            return_total_density = return_total_density,
            return_eta = return_eta,
            store_posterior = store_posterior,
            lambda = lambda,
            regularization = regularization,
            delta = delta,
            sigma = sigma,
            lambda_laplacian = lambda_laplacian,
            lambda_laplacian_effective = lambda_laplacian_effective,
            domain_expansion_steps = domain_expansion_steps,
            domain_expansion_axes = domain_expansion_axes,
            maxit = maxit,
            reltol = reltol,
            x = x,
            y = y,
            z = z,
            optim = optim
        )
    )

    if (return_density) {
        out$density = posterior * as.numeric(pred$density)
        colnames(out$density) = cell_types
    }
    if (return_total_density) {
        out$total_density = as.numeric(pred$density)
    }
    if (return_eta) {
        out$eta = as.numeric(pred$eta)
    }
    if (store_posterior) {
        out$posterior = posterior
    }

    out
}

#' Reconstruct outputs from a spatial linear-density fit
#'
#' Computes derived transcript-level outputs from a fitted
#' `spatial_basis_linear_density_field()` object. This is useful for compact
#' fits where `total_density`, `eta`, or cell-type `density` were not returned
#' during fitting.
#'
#' @inheritParams predict_spatial_basis_log_density_field
#' @param fit Result from `spatial_basis_linear_density_field()`.
#'
#' @return A list containing the requested outputs.
#'
#' @export
predict_spatial_basis_linear_density_field <- function(
    fit,
    transcripts_df = NULL,
    posterior = NULL,
    what = c("total_density", "eta", "density"),
    x = fit$parameters$x %||% "x_location",
    y = fit$parameters$y %||% "y_location",
    z = fit$parameters$z %||% "z_location",
    n_threads = NULL
) {
    what = match.arg(what, several.ok = TRUE)
    if (is.null(fit$eta_basis) || is.null(fit$basis_lattice)) {
        stop("fit must contain eta_basis and basis_lattice.", call. = FALSE)
    }
    if (is.null(n_threads)) {
        n_threads = 0L
    } else if (length(n_threads) != 1L || !is.finite(n_threads) || n_threads < 1) {
        stop("n_threads must be NULL or a positive integer.", call. = FALSE)
    } else {
        n_threads = as.integer(n_threads)
    }

    if (is.null(transcripts_df)) {
        if (is.null(fit$transcript_basis_id) || is.null(fit$transcript_basis_weight)) {
            stop(
                "transcripts_df must be supplied when fit does not contain transcript_basis_id and transcript_basis_weight.",
                call. = FALSE
            )
        }
        basis_id = as.matrix(fit$transcript_basis_id) - 1L
        basis_weight = as.matrix(fit$transcript_basis_weight)
        n = nrow(basis_id)
    } else {
        basis = normalize_spatial_basis(fit$parameters$basis)
        coord_cols = if (basis == "2d") c(x, y) else c(x, y, z)
        missing_cols = setdiff(coord_cols, colnames(transcripts_df))
        if (length(missing_cols) > 0L) {
            stop("Missing coordinate column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
        }
        basis_n_threads = if (n_threads == 0L) NULL else n_threads
        coords = as.matrix(transcripts_df[, coord_cols, drop = FALSE])
        storage.mode(coords) = "double"
        bary = if (basis == "2d") {
            tri_barycentric(coords, s = fit$parameters$s, origin = fit$parameters$origin, n_threads = basis_n_threads)
        } else {
            bcc_barycentric(coords, s = fit$parameters$s, origin = fit$parameters$origin, n_threads = basis_n_threads)
        }
        design = remap_spatial_basis_design(basis_design_from_barycentric(bary), fit$basis_lattice)
        basis_id = design$basis_id - 1L
        basis_weight = design$basis_weight
        n = nrow(basis_id)
    }

    pred = spatial_linear_density_predict_cpp(
        par = fit$eta_basis,
        basis_id = basis_id,
        basis_weight = basis_weight,
        n_threads = n_threads
    )

    out = list()
    if ("total_density" %in% what) {
        out$total_density = as.numeric(pred$density)
    }
    if ("eta" %in% what) {
        out$eta = as.numeric(pred$eta)
    }
    if ("density" %in% what) {
        if (is.null(posterior)) {
            posterior = fit$posterior
        }
        if (is.null(posterior)) {
            stop("posterior must be supplied when requesting density.", call. = FALSE)
        }
        posterior = normalize_posterior_matrix(posterior, n)
        out$density = posterior * as.numeric(pred$density)
        colnames(out$density) = colnames(posterior)
    }

    out
}
