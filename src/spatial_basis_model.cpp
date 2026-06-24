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
    const int n_free = n_cell_types - 1;

    if (basis_weight.nrow() != n || basis_weight.ncol() != n_active) {
        stop("basis_weight must have the same dimensions as basis_id.");
    }
    if (gene_index.size() != n) {
        stop("gene_index must have length nrow(basis_id).");
    }
    if (log_signature.ncol() != n_cell_types) {
        stop("log_signature must have n_cell_types columns.");
    }
    if (par.size() != n_basis * n_free) {
        stop("par has incompatible length.");
    }
    if (edge_to.size() != n_edges) {
        stop("edge_from and edge_to must have the same length.");
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

            for (int k = 0; k < n_free; ++k) {
                f[k] += phi * par[m + n_basis * k];
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

            for (int k = 0; k < n_free; ++k) {
                local_grad[m + n_basis * k] += phi * (p[k] - q[k]);
            }
        }
    }

        objective_by_thread[tid] = local_objective;
    }

    double objective = 0.0;
    NumericVector grad(n_par);

    for (int tid = 0; tid < actual_threads; ++tid) {
        objective += objective_by_thread[tid];
        for (int j = 0; j < n_par; ++j) {
            grad[j] += grad_by_thread[tid][j];
        }
    }

    if (purity_lambda > 0.0 && purity > 0) {
        std::vector<double> logits(n_cell_types);
        std::vector<double> pi(n_cell_types);
        std::vector<double> d_penalty_d_pi(n_cell_types);

        for (int m = 0; m < n_basis; ++m) {
            for (int k = 0; k < n_free; ++k) {
                logits[k] = par[m + n_basis * k];
            }
            logits[n_free] = 0.0;

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

            objective += purity_lambda * penalty;

            for (int k = 0; k < n_free; ++k) {
                grad[m + n_basis * k] += purity_lambda * pi[k] *
                    (d_penalty_d_pi[k] - expected_derivative);
            }
        }
    }

    if (lambda > 0.0) {
        for (int e = 0; e < n_edges; ++e) {
            const int m1 = edge_from[e];
            const int m2 = edge_to[e];

            if (m1 < 0 || m1 >= n_basis || m2 < 0 || m2 >= n_basis) {
                stop("edge indices must be valid basis indices.");
            }

            for (int k = 0; k < n_free; ++k) {
                const int idx1 = m1 + n_basis * k;
                const int idx2 = m2 + n_basis * k;
                const double diff = par[idx1] - par[idx2];

                double penalty = 0.0;
                double derivative = 0.0;

                if (regularization == 0) {
                    penalty = 0.5 * diff * diff;
                    derivative = diff;
                } else if (regularization == 1) {
                    const double abs_diff = std::abs(diff);

                    if (abs_diff <= delta) {
                        penalty = 0.5 * diff * diff;
                        derivative = diff;
                    } else {
                        penalty = delta * (abs_diff - 0.5 * delta);
                        derivative = delta * ((diff >= 0.0) ? 1.0 : -1.0);
                    }
                } else {
                    const double scaled = diff / sigma;
                    const double attenuation = std::exp(-0.5 * scaled * scaled);

                    penalty = 1.0 - attenuation;
                    derivative = attenuation * diff / (sigma * sigma);
                }

                objective += lambda * penalty;
                grad[idx1] += lambda * derivative;
                grad[idx2] -= lambda * derivative;
            }
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
    const int n_cell_types
) {
    const int n = basis_id.nrow();
    const int n_active = basis_id.ncol();
    const int n_free = n_cell_types - 1;

    if (par.size() != n_basis * n_free) {
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

    NumericMatrix prior(n, n_cell_types);
    NumericMatrix posterior(n, n_cell_types);
    NumericMatrix logits(n, n_cell_types);

    const int actual_threads = resolve_threads(n_threads);

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

            for (int k = 0; k < n_free; ++k) {
                f[k] += phi * par[m + n_basis * k];
            }
        }

        const double log_z_prior = log_sum_exp(f);
        const int g = gene_index[i];

        for (int k = 0; k < n_cell_types; ++k) {
            logits(i, k) = f[k];
            prior(i, k) = std::exp(f[k] - log_z_prior);
            log_post[k] = log_signature(g, k) + f[k];
        }

        const double log_z_post = log_sum_exp(log_post);
        for (int k = 0; k < n_cell_types; ++k) {
            posterior(i, k) = std::exp(log_post[k] - log_z_post);
        }
    }
    }

    return List::create(
        _["prior"] = prior,
        _["posterior"] = posterior,
        _["logits"] = logits
    );
}
