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

static inline int log_density_resolve_threads(const int n_threads) {
#ifdef _OPENMP
    if (n_threads <= 0) {
        return omp_get_max_threads();
    }
    return n_threads;
#else
    return 1;
#endif
}

static inline double evaluate_eta_point(
    const int i,
    const IntegerMatrix& basis_id,
    const NumericMatrix& basis_weight,
    const NumericVector& par
) {
    double eta = 0.0;
    const int n_active = basis_id.ncol();
    for (int a = 0; a < n_active; ++a) {
        eta += basis_weight(i, a) * par[basis_id(i, a)];
    }
    return eta;
}

//' @noRd
// [[Rcpp::export]]
List spatial_log_density_objective_cpp(
    const NumericVector& par,
    const IntegerMatrix& obs_basis_id,
    const NumericMatrix& obs_basis_weight,
    const IntegerMatrix& quad_basis_id,
    const NumericMatrix& quad_basis_weight,
    const NumericVector& quad_weight,
    const IntegerVector& edge_from,
    const IntegerVector& edge_to,
    const NumericVector& edge_distance,
    const double lambda,
    const int regularization,
    const double delta,
    const double sigma,
    const double lambda_laplacian,
    const int n_threads
) {
    const int n_obs = obs_basis_id.nrow();
    const int n_quad = quad_basis_id.nrow();
    const int n_active = obs_basis_id.ncol();
    const int n_basis = par.size();
    const int n_edges = edge_from.size();

    if (obs_basis_weight.nrow() != n_obs || obs_basis_weight.ncol() != n_active) {
        stop("obs_basis_weight dimensions do not match obs_basis_id.");
    }
    if (quad_basis_weight.nrow() != n_quad || quad_basis_weight.ncol() != n_active) {
        stop("quad_basis_weight dimensions do not match quad_basis_id.");
    }
    if (quad_weight.size() != n_quad) {
        stop("quad_weight must have one value per quadrature point.");
    }
    if (edge_to.size() != n_edges || edge_distance.size() != n_edges) {
        stop("edge vectors must have matching lengths.");
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
    if (!R_finite(lambda_laplacian) || lambda_laplacian < 0.0) {
        stop("lambda_laplacian must be a non-negative finite number.");
    }

    const int actual_threads = log_density_resolve_threads(n_threads);
    std::vector<double> objective_by_thread(actual_threads, 0.0);
    std::vector< std::vector<double> > grad_by_thread(
        actual_threads,
        std::vector<double>(n_basis, 0.0)
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
        double local_objective = 0.0;
        std::vector<double>& local_grad = grad_by_thread[tid];

#ifdef _OPENMP
#pragma omp for
#endif
        for (int i = 0; i < n_obs; ++i) {
            const double eta = evaluate_eta_point(i, obs_basis_id, obs_basis_weight, par);
            local_objective -= eta;
            for (int a = 0; a < n_active; ++a) {
                local_grad[obs_basis_id(i, a)] -= obs_basis_weight(i, a);
            }
        }

#ifdef _OPENMP
#pragma omp for
#endif
        for (int i = 0; i < n_quad; ++i) {
            const double eta = evaluate_eta_point(i, quad_basis_id, quad_basis_weight, par);
            const double integral_weight = quad_weight[i] * std::exp(eta);
            local_objective += integral_weight;
            for (int a = 0; a < n_active; ++a) {
                local_grad[quad_basis_id(i, a)] += integral_weight * quad_basis_weight(i, a);
            }
        }

        objective_by_thread[tid] = local_objective;
    }

    double objective = 0.0;
    std::vector<double> grad_full(n_basis, 0.0);
    for (int tid = 0; tid < actual_threads; ++tid) {
        objective += objective_by_thread[tid];
        for (int j = 0; j < n_basis; ++j) {
            grad_full[j] += grad_by_thread[tid][j];
        }
    }

    const double inv_n = 1.0 / static_cast<double>(std::max(1, n_obs));
    objective *= inv_n;
    for (int j = 0; j < n_basis; ++j) {
        grad_full[j] *= inv_n;
    }

    if (lambda > 0.0 && n_edges > 0) {
        const double scale = lambda / static_cast<double>(n_edges);
        for (int e = 0; e < n_edges; ++e) {
            const int a = edge_from[e];
            const int b = edge_to[e];
            const double distance = edge_distance[e];
            const double slope = (par[a] - par[b]) / distance;
            const double abs_slope = std::abs(slope);

            double penalty = 0.0;
            double derivative_wrt_slope = 0.0;
            if (regularization == 0) {
                penalty = 0.5 * slope * slope;
                derivative_wrt_slope = slope;
            } else if (regularization == 1) {
                if (abs_slope <= delta) {
                    penalty = 0.5 * slope * slope;
                    derivative_wrt_slope = slope;
                } else {
                    penalty = delta * (abs_slope - 0.5 * delta);
                    derivative_wrt_slope = delta * ((slope >= 0.0) ? 1.0 : -1.0);
                }
            } else {
                const double scaled = slope / sigma;
                const double denom = 1.0 + scaled * scaled;
                const double attenuation = 1.0 / denom;
                penalty = sigma * sigma * (1.0 - attenuation);
                derivative_wrt_slope = 2.0 * slope / (denom * denom);
            }

            objective += scale * penalty;
            const double derivative_wrt_diff = scale * derivative_wrt_slope / distance;
            grad_full[a] += derivative_wrt_diff;
            grad_full[b] -= derivative_wrt_diff;
        }
    }

    if (lambda_laplacian > 0.0 && n_edges > 0) {
        std::vector<double> neighbor_weight_sum(n_basis, 0.0);
        std::vector<double> neighbor_weight(n_edges, 0.0);

        for (int e = 0; e < n_edges; ++e) {
            const int a = edge_from[e];
            const int b = edge_to[e];
            const double distance = edge_distance[e];
            const double w = 1.0 / (distance * distance);
            neighbor_weight[e] = w;
            neighbor_weight_sum[a] += w;
            neighbor_weight_sum[b] += w;
        }

        std::vector<double> laplacian(n_basis, 0.0);
        for (int a = 0; a < n_basis; ++a) {
            laplacian[a] = par[a];
        }
        for (int e = 0; e < n_edges; ++e) {
            const int a = edge_from[e];
            const int b = edge_to[e];
            const double w = neighbor_weight[e];
            if (neighbor_weight_sum[a] > 0.0) {
                laplacian[a] -= w * par[b] / neighbor_weight_sum[a];
            }
            if (neighbor_weight_sum[b] > 0.0) {
                laplacian[b] -= w * par[a] / neighbor_weight_sum[b];
            }
        }

        const double scale = lambda_laplacian / static_cast<double>(n_basis);
        for (int a = 0; a < n_basis; ++a) {
            objective += 0.5 * scale * laplacian[a] * laplacian[a];
            grad_full[a] += scale * laplacian[a];
        }
        for (int e = 0; e < n_edges; ++e) {
            const int a = edge_from[e];
            const int b = edge_to[e];
            const double w = neighbor_weight[e];
            if (neighbor_weight_sum[a] > 0.0) {
                grad_full[b] -= scale * w * laplacian[a] / neighbor_weight_sum[a];
            }
            if (neighbor_weight_sum[b] > 0.0) {
                grad_full[a] -= scale * w * laplacian[b] / neighbor_weight_sum[b];
            }
        }
    }

    NumericVector grad(n_basis);
    for (int j = 0; j < n_basis; ++j) {
        grad[j] = grad_full[j];
    }

    return List::create(
        _["value"] = objective,
        _["gradient"] = grad
    );
}

//' @noRd
// [[Rcpp::export]]
List spatial_log_density_predict_cpp(
    const NumericVector& par,
    const IntegerMatrix& basis_id,
    const NumericMatrix& basis_weight,
    const int n_threads
) {
    const int n = basis_id.nrow();
    NumericVector eta_out(n);
    NumericVector density_out(n);
    const int actual_threads = log_density_resolve_threads(n_threads);

#ifdef _OPENMP
#pragma omp parallel for num_threads(actual_threads)
#endif
    for (int i = 0; i < n; ++i) {
        const double eta = evaluate_eta_point(i, basis_id, basis_weight, par);
        eta_out[i] = eta;
        density_out[i] = std::exp(eta);
    }

    return List::create(
        _["eta"] = eta_out,
        _["density"] = density_out
    );
}
