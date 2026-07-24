// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#endif

static inline int resolve_threads(const int n_threads) {
#ifdef _OPENMP
    if (n_threads <= 0) {
        return omp_get_max_threads();
    }
    return n_threads;
#else
    return 1;
#endif
}

static inline double log_sum_exp(const std::vector<double>& x) {
    double mx = x[0];
    for (std::size_t i = 1; i < x.size(); ++i) {
        if (x[i] > mx) mx = x[i];
    }

    double total = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        total += std::exp(x[i] - mx);
    }

    return mx + std::log(total);
}

static std::vector<double> centered_basis_weights(
    const NumericVector& par,
    const int n_basis,
    const int n_cell_types
) {
    std::vector<double> weights(n_basis * n_cell_types, 0.0);

    for (int m = 0; m < n_basis; ++m) {
        double mean = 0.0;
        for (int k = 0; k < n_cell_types; ++k) {
            mean += par[m + n_basis * k];
        }
        mean /= static_cast<double>(n_cell_types);

        for (int k = 0; k < n_cell_types; ++k) {
            weights[m + n_basis * k] = par[m + n_basis * k] - mean;
        }
    }

    return weights;
}

// C++ backend for spatial_basis_segmentation().
// [[Rcpp::export]]
List spatial_basis_objective_cpp(
    const NumericVector& par,
    const IntegerMatrix& basis_id,
    const NumericMatrix& basis_weight,
    const IntegerVector& gene_index,
    const NumericMatrix& log_signature,
    const IntegerVector& edge_from,
    const IntegerVector& edge_to,
    const NumericVector& edge_distance,
    const double lambda,
    const int regularization,
    const double delta,
    const double sigma,
    const int purity,
    const double purity_lambda,
    const int n_threads,
    const int n_basis,
    const int n_cell_types
) {
    const int n = basis_id.nrow();
    const int n_active = basis_id.ncol();
    const int n_edges = edge_from.size();

    if (basis_weight.nrow() != n || basis_weight.ncol() != n_active) {
        stop("basis_weight must have the same dimensions as basis_id.");
    }
    if (gene_index.size() != n) {
        stop("gene_index must have length nrow(basis_id).");
    }
    if (log_signature.ncol() != n_cell_types) {
        stop("log_signature must have n_cell_types columns.");
    }
    if (par.size() != n_basis * n_cell_types) {
        stop("par has incompatible length.");
    }
    if (edge_to.size() != n_edges) {
        stop("edge_from and edge_to must have the same length.");
    }
    if (edge_distance.size() != n_edges) {
        stop("edge_distance must have the same length as edge_from and edge_to.");
    }
    if (regularization < 0 || regularization > 2) {
        stop("regularization must be 0, 1, or 2.");
    }
    if (!R_finite(delta) || delta <= 0.0) {
        stop("delta must be a positive finite number.");
    }
    if (!R_finite(sigma) || sigma <= 0.0) {
        stop("sigma must be a positive finite number.");
    }
    if (n_threads < 0) {
        stop("n_threads must be NULL or a positive integer.");
    }
    if (purity < 0 || purity > 2) {
        stop("purity must be 0, 1, or 2.");
    }
    if (!R_finite(purity_lambda) || purity_lambda < 0.0) {
        stop("purity_lambda must be a non-negative finite number.");
    }

    for (int i = 0; i < n; ++i) {
        const int g = gene_index[i];
        if (g < 0 || g >= log_signature.nrow()) {
            stop("gene_index contains an out-of-range gene index.");
        }

        for (int a = 0; a < n_active; ++a) {
            const int m = basis_id(i, a);
            if (m < 0 || m >= n_basis) {
                stop("basis_id contains an out-of-range basis index.");
            }
        }
    }

    const int actual_threads = resolve_threads(n_threads);
    const int n_par = par.size();
    const std::vector<double> centered_weights = centered_basis_weights(
        par,
        n_basis,
        n_cell_types
    );
    std::vector<double> objective_by_thread(actual_threads, 0.0);
    std::vector< std::vector<double> > grad_by_thread(
        actual_threads,
        std::vector<double>(n_par, 0.0)
    );

#ifdef _OPENMP
#pragma omp parallel num_threads(actual_threads)
#endif
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        std::vector<double> f(n_cell_types);
        std::vector<double> log_post(n_cell_types);
        std::vector<double> p(n_cell_types);
        std::vector<double> q(n_cell_types);
        double local_objective = 0.0;
        std::vector<double>& local_grad = grad_by_thread[tid];

#ifdef _OPENMP
#pragma omp for
#endif
    for (int i = 0; i < n; ++i) {
        std::fill(f.begin(), f.end(), 0.0);

        for (int a = 0; a < n_active; ++a) {
            const int m = basis_id(i, a);
            const double phi = basis_weight(i, a);

            for (int k = 0; k < n_cell_types; ++k) {
                f[k] += phi * centered_weights[m + n_basis * k];
            }
        }

        const double log_z_prior = log_sum_exp(f);
        for (int k = 0; k < n_cell_types; ++k) {
            p[k] = std::exp(f[k] - log_z_prior);
        }

        const int g = gene_index[i];

        for (int k = 0; k < n_cell_types; ++k) {
            log_post[k] = log_signature(g, k) + f[k];
        }

        const double log_z_post = log_sum_exp(log_post);
        local_objective -= log_z_post - log_z_prior;

        for (int k = 0; k < n_cell_types; ++k) {
            q[k] = std::exp(log_post[k] - log_z_post);
        }

        for (int a = 0; a < n_active; ++a) {
            const int m = basis_id(i, a);
            const double phi = basis_weight(i, a);

            for (int k = 0; k < n_cell_types; ++k) {
                local_grad[m + n_basis * k] += phi * (p[k] - q[k]);
            }
        }
    }

        objective_by_thread[tid] = local_objective;
    }

    double objective = 0.0;
    std::vector<double> grad_full(n_par, 0.0);

    for (int tid = 0; tid < actual_threads; ++tid) {
        objective += objective_by_thread[tid];
        for (int j = 0; j < n_par; ++j) {
            grad_full[j] += grad_by_thread[tid][j];
        }
    }

    const double inv_n = 1.0 / static_cast<double>(n);
    objective *= inv_n;
    for (int j = 0; j < n_par; ++j) {
        grad_full[j] *= inv_n;
    }

    if (purity_lambda > 0.0 && purity > 0) {
        const double purity_scale = purity_lambda / static_cast<double>(n_basis);
        std::vector<double> logits(n_cell_types);
        std::vector<double> pi(n_cell_types);
        std::vector<double> d_penalty_d_pi(n_cell_types);

        for (int m = 0; m < n_basis; ++m) {
            for (int k = 0; k < n_cell_types; ++k) {
                logits[k] = centered_weights[m + n_basis * k];
            }

            const double log_z = log_sum_exp(logits);
            for (int k = 0; k < n_cell_types; ++k) {
                pi[k] = std::exp(logits[k] - log_z);
            }

            double penalty = 0.0;

            if (purity == 1) {
                for (int k = 0; k < n_cell_types; ++k) {
                    penalty -= pi[k] * std::log(pi[k]);
                    d_penalty_d_pi[k] = -(std::log(pi[k]) + 1.0);
                }
            } else {
                double sum_pi2 = 0.0;
                for (int k = 0; k < n_cell_types; ++k) {
                    sum_pi2 += pi[k] * pi[k];
                    d_penalty_d_pi[k] = -2.0 * pi[k];
                }
                penalty = 1.0 - sum_pi2;
            }

            double expected_derivative = 0.0;
            for (int k = 0; k < n_cell_types; ++k) {
                expected_derivative += pi[k] * d_penalty_d_pi[k];
            }

            objective += purity_scale * penalty;

            for (int k = 0; k < n_cell_types; ++k) {
                grad_full[m + n_basis * k] += purity_scale * pi[k] *
                    (d_penalty_d_pi[k] - expected_derivative);
            }
        }
    }

    if (lambda > 0.0 && n_edges > 0) {
        const double edge_scale = lambda / static_cast<double>(n_edges);
        for (int e = 0; e < n_edges; ++e) {
            const int m1 = edge_from[e];
            const int m2 = edge_to[e];
            const double distance = edge_distance[e];

            if (m1 < 0 || m1 >= n_basis || m2 < 0 || m2 >= n_basis) {
                stop("edge indices must be valid basis indices.");
            }
            if (!R_finite(distance) || distance <= 0.0) {
                stop("edge_distance must contain positive finite values.");
            }

            for (int k = 0; k < n_cell_types; ++k) {
                const int idx1 = m1 + n_basis * k;
                const int idx2 = m2 + n_basis * k;
                const double diff = centered_weights[idx1] - centered_weights[idx2];
                const double slope = diff / distance;

                double penalty = 0.0;
                double derivative_wrt_slope = 0.0;

                if (regularization == 0) {
                    penalty = 0.5 * slope * slope;
                    derivative_wrt_slope = slope;
                } else if (regularization == 1) {
                    const double abs_slope = std::abs(slope);

                    if (abs_slope <= delta) {
                        penalty = 0.5 * slope * slope;
                        derivative_wrt_slope = slope;
                    } else {
                        penalty = delta * (abs_slope - 0.5 * delta);
                        derivative_wrt_slope = delta * ((slope >= 0.0) ? 1.0 : -1.0);
                    }
                } else {
                    const double scaled = slope / sigma;
                    const double attenuation = std::exp(-0.5 * scaled * scaled);

                    penalty = sigma * sigma * (1.0 - attenuation);
                    derivative_wrt_slope = attenuation * slope;
                }

                const double derivative = derivative_wrt_slope / distance;
                objective += edge_scale * penalty;
                grad_full[idx1] += edge_scale * derivative;
                grad_full[idx2] -= edge_scale * derivative;
            }
        }
    }

    NumericVector grad(n_par);
    for (int m = 0; m < n_basis; ++m) {
        double mean_grad = 0.0;
        for (int k = 0; k < n_cell_types; ++k) {
            mean_grad += grad_full[m + n_basis * k];
        }
        mean_grad /= static_cast<double>(n_cell_types);

        for (int k = 0; k < n_cell_types; ++k) {
            grad[m + n_basis * k] = grad_full[m + n_basis * k] - mean_grad;
        }
    }

    return List::create(
        _["value"] = objective,
        _["gradient"] = grad
    );
}

// C++ backend for spatial_basis_segmentation().
// [[Rcpp::export]]
List spatial_basis_predict_cpp(
    const NumericVector& par,
    const IntegerMatrix& basis_id,
    const NumericMatrix& basis_weight,
    const IntegerVector& gene_index,
    const NumericMatrix& log_signature,
    const int n_basis,
    const int n_threads,
    const int n_cell_types,
    const bool return_prior,
    const bool return_posterior,
    const bool return_logits,
    const bool return_labels,
    const bool return_max_posterior
) {
    const int n = basis_id.nrow();
    const int n_active = basis_id.ncol();

    if (par.size() != n_basis * n_cell_types) {
        stop("par has incompatible length.");
    }
    if (n_threads < 0) {
        stop("n_threads must be NULL or a positive integer.");
    }

    for (int i = 0; i < n; ++i) {
        const int g = gene_index[i];
        if (g < 0 || g >= log_signature.nrow()) {
            stop("gene_index contains an out-of-range gene index.");
        }

        for (int a = 0; a < n_active; ++a) {
            const int m = basis_id(i, a);
            if (m < 0 || m >= n_basis) {
                stop("basis_id contains an out-of-range basis index.");
            }
        }
    }

    NumericMatrix prior;
    NumericMatrix posterior;
    NumericMatrix logits;
    IntegerVector labels;
    NumericVector max_posterior;

    if (return_prior) {
        prior = NumericMatrix(n, n_cell_types);
    }
    if (return_posterior) {
        posterior = NumericMatrix(n, n_cell_types);
    }
    if (return_logits) {
        logits = NumericMatrix(n, n_cell_types);
    }
    if (return_labels) {
        labels = IntegerVector(n);
    }
    if (return_max_posterior) {
        max_posterior = NumericVector(n);
    }

    const int actual_threads = resolve_threads(n_threads);
    const std::vector<double> centered_weights = centered_basis_weights(
        par,
        n_basis,
        n_cell_types
    );

#ifdef _OPENMP
#pragma omp parallel num_threads(actual_threads)
#endif
    {
        std::vector<double> f(n_cell_types);
        std::vector<double> log_post(n_cell_types);

#ifdef _OPENMP
#pragma omp for
#endif
    for (int i = 0; i < n; ++i) {
        std::fill(f.begin(), f.end(), 0.0);

        for (int a = 0; a < n_active; ++a) {
            const int m = basis_id(i, a);
            const double phi = basis_weight(i, a);

            for (int k = 0; k < n_cell_types; ++k) {
                f[k] += phi * centered_weights[m + n_basis * k];
            }
        }

        const double log_z_prior = log_sum_exp(f);
        const int g = gene_index[i];

        for (int k = 0; k < n_cell_types; ++k) {
            if (return_logits) {
                logits(i, k) = f[k];
            }
            if (return_prior) {
                prior(i, k) = std::exp(f[k] - log_z_prior);
            }
            log_post[k] = log_signature(g, k) + f[k];
        }

        const double log_z_post = log_sum_exp(log_post);
        int best_k = 0;
        double best_q = R_NegInf;

        for (int k = 0; k < n_cell_types; ++k) {
            const double q = std::exp(log_post[k] - log_z_post);
            if (return_posterior) {
                posterior(i, k) = q;
            }
            if (q > best_q) {
                best_q = q;
                best_k = k;
            }
        }

        if (return_labels) {
            labels[i] = best_k + 1;
        }
        if (return_max_posterior) {
            max_posterior[i] = best_q;
        }
    }
    }

    List out;
    if (return_prior) {
        out["prior"] = prior;
    }
    if (return_posterior) {
        out["posterior"] = posterior;
    }
    if (return_logits) {
        out["logits"] = logits;
    }
    if (return_labels) {
        out["labels"] = labels;
    }
    if (return_max_posterior) {
        out["max_posterior"] = max_posterior;
    }

    return out;
}

// C++ backend for streamed signature refinement.
// [[Rcpp::export]]
List spatial_basis_posterior_counts_cpp(
    const NumericVector& par,
    const IntegerMatrix& basis_id,
    const NumericMatrix& basis_weight,
    const IntegerVector& gene_index,
    const NumericMatrix& log_signature,
    const int n_basis,
    const int n_threads,
    const int n_cell_types,
    const double min_posterior
) {
    const int n = basis_id.nrow();
    const int n_active = basis_id.ncol();
    const int n_genes = log_signature.nrow();

    if (par.size() != n_basis * n_cell_types) {
        stop("par has incompatible length.");
    }
    if (n_threads < 0) {
        stop("n_threads must be NULL or a positive integer.");
    }
    if (!R_finite(min_posterior) || min_posterior < 0.0 || min_posterior > 1.0) {
        stop("min_posterior must be a finite number in [0, 1].");
    }

    for (int i = 0; i < n; ++i) {
        const int g = gene_index[i];
        if (g < 0 || g >= n_genes) {
            stop("gene_index contains an out-of-range gene index.");
        }

        for (int a = 0; a < n_active; ++a) {
            const int m = basis_id(i, a);
            if (m < 0 || m >= n_basis) {
                stop("basis_id contains an out-of-range basis index.");
            }
        }
    }

    const int actual_threads = resolve_threads(n_threads);
    const std::vector<double> centered_weights = centered_basis_weights(
        par,
        n_basis,
        n_cell_types
    );
    std::vector< std::vector<double> > counts_by_thread(
        actual_threads,
        std::vector<double>(n_genes * n_cell_types, 0.0)
    );

#ifdef _OPENMP
#pragma omp parallel num_threads(actual_threads)
#endif
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        std::vector<double> f(n_cell_types);
        std::vector<double> log_post(n_cell_types);
        std::vector<double>& local_counts = counts_by_thread[tid];

#ifdef _OPENMP
#pragma omp for
#endif
    for (int i = 0; i < n; ++i) {
        std::fill(f.begin(), f.end(), 0.0);

        for (int a = 0; a < n_active; ++a) {
            const int m = basis_id(i, a);
            const double phi = basis_weight(i, a);

            for (int k = 0; k < n_cell_types; ++k) {
                f[k] += phi * centered_weights[m + n_basis * k];
            }
        }

        const int g = gene_index[i];

        for (int k = 0; k < n_cell_types; ++k) {
            log_post[k] = log_signature(g, k) + f[k];
        }

        const double log_z_post = log_sum_exp(log_post);
        for (int k = 0; k < n_cell_types; ++k) {
            const double q = std::exp(log_post[k] - log_z_post);
            if (q >= min_posterior) {
                local_counts[g + n_genes * k] += q;
            }
        }
    }
    }

    NumericMatrix counts(n_genes, n_cell_types);
    NumericVector effective_counts(n_cell_types);

    for (int tid = 0; tid < actual_threads; ++tid) {
        const std::vector<double>& local_counts = counts_by_thread[tid];
        for (int k = 0; k < n_cell_types; ++k) {
            for (int g = 0; g < n_genes; ++g) {
                const double value = local_counts[g + n_genes * k];
                counts(g, k) += value;
                effective_counts[k] += value;
            }
        }
    }

    return List::create(
        _["counts"] = counts,
        _["effective_counts"] = effective_counts
    );
}

// C++ backend for density-basis posterior support without storing marginals.
// [[Rcpp::export]]
NumericMatrix spatial_basis_density_support_cpp(
    const NumericVector& par,
    const IntegerMatrix& segmentation_basis_id,
    const NumericMatrix& segmentation_basis_weight,
    const IntegerVector& gene_index,
    const NumericMatrix& log_signature,
    const IntegerMatrix& density_basis_id,
    const NumericMatrix& density_basis_weight,
    const int n_segmentation_basis,
    const int n_density_basis,
    const int n_threads,
    const int n_cell_types,
    const bool hard_max
) {
    const int n = segmentation_basis_id.nrow();
    const int n_active_segmentation = segmentation_basis_id.ncol();
    const int n_active_density = density_basis_id.ncol();

    if (segmentation_basis_weight.nrow() != n || segmentation_basis_weight.ncol() != n_active_segmentation) {
        stop("segmentation_basis_weight must have the same dimensions as segmentation_basis_id.");
    }
    if (density_basis_weight.nrow() != n || density_basis_weight.ncol() != n_active_density) {
        stop("density_basis_weight must have the same dimensions as density_basis_id.");
    }
    if (gene_index.size() != n) {
        stop("gene_index must have length nrow(segmentation_basis_id).");
    }
    if (log_signature.ncol() != n_cell_types) {
        stop("log_signature must have n_cell_types columns.");
    }
    if (par.size() != n_segmentation_basis * n_cell_types) {
        stop("par has incompatible length.");
    }
    if (n_threads < 0) {
        stop("n_threads must be NULL or a positive integer.");
    }

    for (int i = 0; i < n; ++i) {
        const int g = gene_index[i];
        if (g < 0 || g >= log_signature.nrow()) {
            stop("gene_index contains an out-of-range gene index.");
        }

        for (int a = 0; a < n_active_segmentation; ++a) {
            const int m = segmentation_basis_id(i, a);
            if (m < 0 || m >= n_segmentation_basis) {
                stop("segmentation_basis_id contains an out-of-range basis index.");
            }
        }

        for (int a = 0; a < n_active_density; ++a) {
            const int m = density_basis_id(i, a);
            if (m != NA_INTEGER && (m < 0 || m >= n_density_basis)) {
                stop("density_basis_id contains an out-of-range basis index.");
            }
        }
    }

    const int actual_threads = resolve_threads(n_threads);
    const std::vector<double> centered_weights = centered_basis_weights(
        par,
        n_segmentation_basis,
        n_cell_types
    );
    std::vector< std::vector<double> > support_by_thread(
        actual_threads,
        std::vector<double>(n_density_basis * n_cell_types, 0.0)
    );

#ifdef _OPENMP
#pragma omp parallel num_threads(actual_threads)
#endif
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        std::vector<double> f(n_cell_types);
        std::vector<double> log_post(n_cell_types);
        std::vector<double> posterior(n_cell_types);
        std::vector<double>& local_support = support_by_thread[tid];

#ifdef _OPENMP
#pragma omp for
#endif
    for (int i = 0; i < n; ++i) {
        std::fill(f.begin(), f.end(), 0.0);

        for (int a = 0; a < n_active_segmentation; ++a) {
            const int m = segmentation_basis_id(i, a);
            const double phi = segmentation_basis_weight(i, a);

            for (int k = 0; k < n_cell_types; ++k) {
                f[k] += phi * centered_weights[m + n_segmentation_basis * k];
            }
        }

        const int g = gene_index[i];
        for (int k = 0; k < n_cell_types; ++k) {
            log_post[k] = log_signature(g, k) + f[k];
        }

        const double log_z_post = log_sum_exp(log_post);
        int best_k = 0;
        double best_log_post = log_post[0];
        for (int k = 0; k < n_cell_types; ++k) {
            if (log_post[k] > best_log_post) {
                best_log_post = log_post[k];
                best_k = k;
            }
            posterior[k] = std::exp(log_post[k] - log_z_post);
        }

        for (int a = 0; a < n_active_density; ++a) {
            const int m = density_basis_id(i, a);
            if (m == NA_INTEGER) {
                continue;
            }
            const double phi = density_basis_weight(i, a);
            if (!R_finite(phi) || phi == 0.0) {
                continue;
            }

            if (hard_max) {
                local_support[m + n_density_basis * best_k] += phi;
            } else {
                for (int k = 0; k < n_cell_types; ++k) {
                    local_support[m + n_density_basis * k] += phi * posterior[k];
                }
            }
        }
    }
    }

    NumericMatrix support(n_density_basis, n_cell_types);
    for (int tid = 0; tid < actual_threads; ++tid) {
        const std::vector<double>& local_support = support_by_thread[tid];
        for (int k = 0; k < n_cell_types; ++k) {
            for (int m = 0; m < n_density_basis; ++m) {
                support(m, k) += local_support[m + n_density_basis * k];
            }
        }
    }

    return support;
}

// Pairwise overlap summaries for spatial prior fields.
// [[Rcpp::export]]
List spatial_basis_signature_overlap_cpp(
    const NumericMatrix& spatial_prior,
    const NumericVector& weights,
    const int n_threads
) {
    const int n = spatial_prior.nrow();
    const int K = spatial_prior.ncol();
    if (weights.size() != n) {
        stop("weights must have one entry per row of spatial_prior.");
    }
    if (K < 2) {
        stop("spatial_prior must contain at least two columns.");
    }
    if (n_threads < 0) {
        stop("n_threads must be NULL or a positive integer.");
    }

    double total_weight = 0.0;
    for (int i = 0; i < n; ++i) {
        const double wi = weights[i];
        if (!R_finite(wi) || wi < 0.0) {
            stop("weights must contain only non-negative finite values.");
        }
        total_weight += wi;
        for (int k = 0; k < K; ++k) {
            if (!R_finite(spatial_prior(i, k))) {
                stop("spatial_prior must contain only finite values.");
            }
        }
    }
    if (total_weight <= 0.0) {
        stop("weights must have positive total mass.");
    }

    const int P = K * (K - 1) / 2;
    IntegerVector pair_a(P);
    IntegerVector pair_b(P);
    int p = 0;
    for (int a = 0; a < K - 1; ++a) {
        for (int b = a + 1; b < K; ++b) {
            pair_a[p] = a + 1;
            pair_b[p] = b + 1;
            ++p;
        }
    }

    NumericVector mass(K);
    NumericVector norm2(K);
    NumericVector shared(P);
    NumericVector cross(P);
    const int actual_threads = resolve_threads(n_threads);

#ifdef _OPENMP
#pragma omp parallel num_threads(actual_threads)
    {
        std::vector<double> local_mass(K, 0.0);
        std::vector<double> local_norm2(K, 0.0);
        std::vector<double> local_shared(P, 0.0);
        std::vector<double> local_cross(P, 0.0);

#pragma omp for schedule(static)
        for (int i = 0; i < n; ++i) {
            const double wi = weights[i];
            for (int k = 0; k < K; ++k) {
                const double pk = spatial_prior(i, k);
                local_mass[k] += wi * pk;
                local_norm2[k] += wi * pk * pk;
            }
            int q = 0;
            for (int a = 0; a < K - 1; ++a) {
                const double pa = spatial_prior(i, a);
                for (int b = a + 1; b < K; ++b) {
                    const double pb = spatial_prior(i, b);
                    local_shared[q] += wi * std::min(pa, pb);
                    local_cross[q] += wi * pa * pb;
                    ++q;
                }
            }
        }

#pragma omp critical
        {
            for (int k = 0; k < K; ++k) {
                mass[k] += local_mass[k];
                norm2[k] += local_norm2[k];
            }
            for (int q = 0; q < P; ++q) {
                shared[q] += local_shared[q];
                cross[q] += local_cross[q];
            }
        }
    }
#else
    for (int i = 0; i < n; ++i) {
        const double wi = weights[i];
        for (int k = 0; k < K; ++k) {
            const double pk = spatial_prior(i, k);
            mass[k] += wi * pk;
            norm2[k] += wi * pk * pk;
        }
        int q = 0;
        for (int a = 0; a < K - 1; ++a) {
            const double pa = spatial_prior(i, a);
            for (int b = a + 1; b < K; ++b) {
                const double pb = spatial_prior(i, b);
                shared[q] += wi * std::min(pa, pb);
                cross[q] += wi * pa * pb;
                ++q;
            }
        }
    }
#endif

    return List::create(
        _["pair_a"] = pair_a,
        _["pair_b"] = pair_b,
        _["mass"] = mass,
        _["norm2"] = norm2,
        _["shared"] = shared,
        _["cross"] = cross
    );
}
