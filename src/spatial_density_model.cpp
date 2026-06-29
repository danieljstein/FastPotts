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

static inline int density_resolve_threads(const int n_threads) {
#ifdef _OPENMP
    if (n_threads <= 0) {
        return omp_get_max_threads();
    }
    return n_threads;
#else
    return 1;
#endif
}

static inline double sigmoid_neg_stable(const double g) {
    if (g >= 0.0) {
        const double z = std::exp(-g);
        return z / (1.0 + z);
    }
    const double z = std::exp(g);
    return 1.0 / (1.0 + z);
}

static inline double attenuation_and_h(
    const double g,
    const double floor_value,
    double& h
) {
    const double sig_neg = sigmoid_neg_stable(g);
    const double attenuation = floor_value + (1.0 - floor_value) * sig_neg;
    const double d_attenuation_d_g = -(1.0 - floor_value) * sig_neg * (1.0 - sig_neg);
    h = d_attenuation_d_g / attenuation;
    return attenuation;
}

static inline void evaluate_density_point(
    const int i,
    const int k,
    const IntegerMatrix& basis_id,
    const NumericMatrix& basis_weight,
    const IntegerMatrix& edge_id,
    const IntegerVector& pair_from,
    const IntegerVector& pair_to,
    const NumericVector& par,
    const int n_basis,
    const int n_edge_coef,
    const int n_cell_types,
    const NumericVector& density_floor,
    double& eta,
    double& g,
    double& attenuation,
    double& h
) {
    const int n_active = basis_id.ncol();
    const int n_pairs = edge_id.ncol();
    const int eta_offset = 0;
    const int node_offset = n_basis * n_cell_types;
    const int edge_offset = 2 * n_basis * n_cell_types;

    eta = 0.0;
    g = 0.0;

    for (int a = 0; a < n_active; ++a) {
        const int m = basis_id(i, a);
        const double phi = basis_weight(i, a);
        eta += phi * par[eta_offset + m + n_basis * k];
        g += phi * par[node_offset + m + n_basis * k];
    }

    for (int p = 0; p < n_pairs; ++p) {
        const int e = edge_id(i, p);
        const double feature = basis_weight(i, pair_from[p]) *
            basis_weight(i, pair_to[p]);
        g += feature * par[edge_offset + e + n_edge_coef * k];
    }

    attenuation = attenuation_and_h(g, density_floor[k], h);
}

//' @noRd
// [[Rcpp::export]]
List spatial_density_objective_cpp(
    const NumericVector& par,
    const IntegerMatrix& obs_basis_id,
    const NumericMatrix& obs_basis_weight,
    const IntegerMatrix& obs_edge_id,
    const NumericMatrix& obs_weight,
    const IntegerMatrix& quad_basis_id,
    const NumericMatrix& quad_basis_weight,
    const IntegerMatrix& quad_edge_id,
    const NumericVector& quad_weight,
    const IntegerVector& pair_from,
    const IntegerVector& pair_to,
    const IntegerVector& smooth_edge_from,
    const IntegerVector& smooth_edge_to,
    const NumericVector& smooth_edge_distance,
    const NumericVector& density_floor,
    const double lambda_eta,
    const double lambda_s_smooth,
    const double lambda_s_prior,
    const double s_prior_mean,
    const double lambda_edge,
    const int n_basis,
    const int n_edge_coef,
    const int n_cell_types,
    const int n_threads
) {
    const int n_obs = obs_basis_id.nrow();
    const int n_quad = quad_basis_id.nrow();
    const int n_active = obs_basis_id.ncol();
    const int n_pairs = obs_edge_id.ncol();
    const int n_smooth_edges = smooth_edge_from.size();
    const int n_par_expected = (2 * n_basis + n_edge_coef) * n_cell_types;

    if (par.size() != n_par_expected) {
        stop("par has incompatible length.");
    }
    if (obs_basis_weight.nrow() != n_obs || obs_basis_weight.ncol() != n_active) {
        stop("obs_basis_weight dimensions do not match obs_basis_id.");
    }
    if (quad_basis_weight.nrow() != n_quad || quad_basis_weight.ncol() != n_active) {
        stop("quad_basis_weight dimensions do not match quad_basis_id.");
    }
    if (obs_edge_id.nrow() != n_obs || quad_edge_id.nrow() != n_quad ||
        quad_edge_id.ncol() != n_pairs) {
        stop("edge_id dimensions are incompatible.");
    }
    if (obs_weight.nrow() != n_obs || obs_weight.ncol() != n_cell_types) {
        stop("obs_weight must have n_obs rows and n_cell_types columns.");
    }
    if (quad_weight.size() != n_quad) {
        stop("quad_weight must have one value per quadrature point.");
    }
    if (density_floor.size() != n_cell_types) {
        stop("density_floor must have one value per cell type.");
    }

    for (int k = 0; k < n_cell_types; ++k) {
        if (!R_finite(density_floor[k]) || density_floor[k] < 0.0 || density_floor[k] >= 1.0) {
            stop("density_floor values must be in [0, 1).");
        }
    }

    const int eta_offset = 0;
    const int node_offset = n_basis * n_cell_types;
    const int edge_offset = 2 * n_basis * n_cell_types;
    const int actual_threads = density_resolve_threads(n_threads);
    std::vector<double> objective_by_thread(actual_threads, 0.0);
    std::vector< std::vector<double> > grad_by_thread(
        actual_threads,
        std::vector<double>(par.size(), 0.0)
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
            for (int k = 0; k < n_cell_types; ++k) {
                const double weight = obs_weight(i, k);
                if (weight == 0.0) {
                    continue;
                }

                double eta, g, attenuation, h;
                evaluate_density_point(
                    i, k, obs_basis_id, obs_basis_weight, obs_edge_id,
                    pair_from, pair_to, par, n_basis, n_edge_coef,
                    n_cell_types, density_floor, eta, g, attenuation, h
                );

                local_objective -= weight * (eta + std::log(attenuation));

                for (int a = 0; a < n_active; ++a) {
                    const int m = obs_basis_id(i, a);
                    const double phi = obs_basis_weight(i, a);
                    local_grad[eta_offset + m + n_basis * k] -= weight * phi;
                    local_grad[node_offset + m + n_basis * k] -= weight * h * phi;
                }

                for (int p = 0; p < n_pairs; ++p) {
                    const int e = obs_edge_id(i, p);
                    const double feature = obs_basis_weight(i, pair_from[p]) *
                        obs_basis_weight(i, pair_to[p]);
                    local_grad[edge_offset + e + n_edge_coef * k] -= weight * h * feature;
                }
            }
        }

#ifdef _OPENMP
#pragma omp for
#endif
        for (int i = 0; i < n_quad; ++i) {
            const double q_weight = quad_weight[i];
            for (int k = 0; k < n_cell_types; ++k) {
                double eta, g, attenuation, h;
                evaluate_density_point(
                    i, k, quad_basis_id, quad_basis_weight, quad_edge_id,
                    pair_from, pair_to, par, n_basis, n_edge_coef,
                    n_cell_types, density_floor, eta, g, attenuation, h
                );

                const double rho = std::exp(eta) * attenuation;
                const double integral_weight = q_weight * rho;
                local_objective += integral_weight;

                for (int a = 0; a < n_active; ++a) {
                    const int m = quad_basis_id(i, a);
                    const double phi = quad_basis_weight(i, a);
                    local_grad[eta_offset + m + n_basis * k] += integral_weight * phi;
                    local_grad[node_offset + m + n_basis * k] += integral_weight * h * phi;
                }

                for (int p = 0; p < n_pairs; ++p) {
                    const int e = quad_edge_id(i, p);
                    const double feature = quad_basis_weight(i, pair_from[p]) *
                        quad_basis_weight(i, pair_to[p]);
                    local_grad[edge_offset + e + n_edge_coef * k] += integral_weight * h * feature;
                }
            }
        }

        objective_by_thread[tid] = local_objective;
    }

    double objective = 0.0;
    std::vector<double> grad_full(par.size(), 0.0);
    for (int tid = 0; tid < actual_threads; ++tid) {
        objective += objective_by_thread[tid];
        for (R_xlen_t j = 0; j < par.size(); ++j) {
            grad_full[j] += grad_by_thread[tid][j];
        }
    }

    double total_weight = 0.0;
    for (int i = 0; i < n_obs; ++i) {
        for (int k = 0; k < n_cell_types; ++k) {
            total_weight += obs_weight(i, k);
        }
    }
    if (total_weight <= 0.0) {
        stop("obs_weight must contain positive total mass.");
    }

    const double inv_total_weight = 1.0 / total_weight;
    objective *= inv_total_weight;
    for (R_xlen_t j = 0; j < par.size(); ++j) {
        grad_full[j] *= inv_total_weight;
    }

    if ((lambda_eta > 0.0 || lambda_s_smooth > 0.0) && n_smooth_edges > 0) {
        for (int e = 0; e < n_smooth_edges; ++e) {
            const int m1 = smooth_edge_from[e];
            const int m2 = smooth_edge_to[e];
            const double distance = smooth_edge_distance[e];
            const double edge_scale_eta = lambda_eta / static_cast<double>(n_smooth_edges);
            const double edge_scale_s = lambda_s_smooth / static_cast<double>(n_smooth_edges);

            for (int k = 0; k < n_cell_types; ++k) {
                if (lambda_eta > 0.0) {
                    const int idx1 = eta_offset + m1 + n_basis * k;
                    const int idx2 = eta_offset + m2 + n_basis * k;
                    const double slope = (par[idx1] - par[idx2]) / distance;
                    const double derivative = edge_scale_eta * slope / distance;
                    objective += edge_scale_eta * 0.5 * slope * slope;
                    grad_full[idx1] += derivative;
                    grad_full[idx2] -= derivative;
                }

                if (lambda_s_smooth > 0.0) {
                    const int idx1 = node_offset + m1 + n_basis * k;
                    const int idx2 = node_offset + m2 + n_basis * k;
                    const double slope = (par[idx1] - par[idx2]) / distance;
                    const double derivative = edge_scale_s * slope / distance;
                    objective += edge_scale_s * 0.5 * slope * slope;
                    grad_full[idx1] += derivative;
                    grad_full[idx2] -= derivative;
                }
            }
        }
    }

    if (lambda_s_prior > 0.0) {
        const double scale = lambda_s_prior / static_cast<double>(n_basis * n_cell_types);
        for (int k = 0; k < n_cell_types; ++k) {
            for (int m = 0; m < n_basis; ++m) {
                const int idx = node_offset + m + n_basis * k;
                const double diff = par[idx] - s_prior_mean;
                objective += scale * 0.5 * diff * diff;
                grad_full[idx] += scale * diff;
            }
        }
    }

    if (lambda_edge > 0.0) {
        const double scale = lambda_edge / static_cast<double>(n_edge_coef * n_cell_types);
        for (int k = 0; k < n_cell_types; ++k) {
            for (int e = 0; e < n_edge_coef; ++e) {
                const int idx = edge_offset + e + n_edge_coef * k;
                objective += scale * 0.5 * par[idx] * par[idx];
                grad_full[idx] += scale * par[idx];
            }
        }
    }

    NumericVector grad(par.size());
    for (R_xlen_t j = 0; j < par.size(); ++j) {
        grad[j] = grad_full[j];
    }

    return List::create(
        _["value"] = objective,
        _["gradient"] = grad
    );
}

//' @noRd
// [[Rcpp::export]]
List spatial_density_predict_cpp(
    const NumericVector& par,
    const IntegerMatrix& basis_id,
    const NumericMatrix& basis_weight,
    const IntegerMatrix& edge_id,
    const IntegerVector& pair_from,
    const IntegerVector& pair_to,
    const NumericVector& density_floor,
    const int n_basis,
    const int n_edge_coef,
    const int n_cell_types,
    const int n_threads
) {
    const int n = basis_id.nrow();
    NumericMatrix eta_out(n, n_cell_types);
    NumericMatrix g_out(n, n_cell_types);
    NumericMatrix attenuation_out(n, n_cell_types);
    NumericMatrix density_out(n, n_cell_types);

    const int actual_threads = density_resolve_threads(n_threads);
#ifdef _OPENMP
#pragma omp parallel for num_threads(actual_threads)
#endif
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n_cell_types; ++k) {
            double eta, g, attenuation, h;
            evaluate_density_point(
                i, k, basis_id, basis_weight, edge_id, pair_from, pair_to,
                par, n_basis, n_edge_coef, n_cell_types, density_floor,
                eta, g, attenuation, h
            );
            eta_out(i, k) = eta;
            g_out(i, k) = g;
            attenuation_out(i, k) = attenuation;
            density_out(i, k) = std::exp(eta) * attenuation;
        }
    }

    return List::create(
        _["eta"] = eta_out,
        _["g"] = g_out,
        _["attenuation"] = attenuation_out,
        _["density"] = density_out
    );
}
