// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
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

static inline double log_density_factorial(const int n) {
    double out = 1.0;
    for (int i = 2; i <= n; ++i) {
        out *= static_cast<double>(i);
    }
    return out;
}

static inline double simplex_log_density_direct(
    const std::vector<double>& h,
    const double volume,
    std::vector<double>& grad
) {
    const int n = h.size();
    const double shift = *std::max_element(h.begin(), h.end());
    std::vector<double> r(n);
    for (int i = 0; i < n; ++i) {
        r[i] = h[i] - shift;
    }

    std::vector<double> term(n, 0.0);
    double dd = 0.0;
    for (int i = 0; i < n; ++i) {
        double denom = 1.0;
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                denom *= r[i] - r[j];
            }
        }
        term[i] = std::exp(r[i]) / denom;
        dd += term[i];
    }

    const double scale = volume * log_density_factorial(n - 1) * std::exp(shift);
    for (int k = 0; k < n; ++k) {
        double inv_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (j != k) {
                inv_sum += 1.0 / (r[k] - r[j]);
            }
        }
        double ddk = term[k] * (1.0 - inv_sum);
        for (int i = 0; i < n; ++i) {
            if (i != k) {
                ddk += term[i] / (r[i] - r[k]);
            }
        }
        grad[k] = scale * ddk;
    }

    return scale * dd;
}

static void enumerate_taylor_order(
    const std::vector<double>& r,
    const int pos,
    const int remaining,
    std::vector<int>& alpha,
    double monomial,
    double& h_m,
    std::vector<double>& grad_poly
) {
    const int n = r.size();
    if (pos == n - 1) {
        alpha[pos] = remaining;
        double term = monomial;
        if (remaining > 0) {
            term *= std::pow(r[pos], remaining);
        }
        h_m += term;
        for (int i = 0; i < n; ++i) {
            grad_poly[i] += static_cast<double>(alpha[i] + 1) * term;
        }
        return;
    }

    double power = 1.0;
    for (int a = 0; a <= remaining; ++a) {
        alpha[pos] = a;
        enumerate_taylor_order(
            r,
            pos + 1,
            remaining - a,
            alpha,
            monomial * power,
            h_m,
            grad_poly
        );
        power *= r[pos];
    }
}

static inline double simplex_log_density_taylor(
    const std::vector<double>& h,
    const double volume,
    std::vector<double>& grad,
    const double tol,
    const int max_terms
) {
    const int n = h.size();
    double center = 0.0;
    for (int i = 0; i < n; ++i) {
        center += h[i];
    }
    center /= static_cast<double>(n);

    std::vector<double> r(n);
    for (int i = 0; i < n; ++i) {
        r[i] = h[i] - center;
    }

    std::fill(grad.begin(), grad.end(), 0.0);
    double value_scaled = 0.0;
    std::vector<int> alpha(n, 0);
    const double gamma_n = log_density_factorial(n - 1);

    for (int m = 0; m <= max_terms; ++m) {
        double h_m = 0.0;
        std::vector<double> grad_poly(n, 0.0);
        enumerate_taylor_order(r, 0, m, alpha, 1.0, h_m, grad_poly);

        double denom_value = 1.0;
        for (int q = 0; q < m; ++q) {
            denom_value *= static_cast<double>(n + q);
        }
        const double coeff_value = gamma_n / (gamma_n * denom_value);
        const double coeff_grad = coeff_value / static_cast<double>(n + m);

        const double value_term = coeff_value * h_m;
        value_scaled += value_term;

        double max_grad_term = std::abs(value_term);
        for (int i = 0; i < n; ++i) {
            const double grad_term = coeff_grad * grad_poly[i];
            grad[i] += grad_term;
            max_grad_term = std::max(max_grad_term, std::abs(grad_term));
        }

        if (m >= 4 && max_grad_term <= tol * std::max(1.0, std::abs(value_scaled))) {
            break;
        }
    }

    const double scale = volume * std::exp(center);
    for (int i = 0; i < n; ++i) {
        grad[i] *= scale;
    }
    return scale * value_scaled;
}

static void gauss_legendre_unit_interval(
    const int q,
    std::vector<double>& nodes,
    std::vector<double>& weights
) {
    nodes.assign(q, 0.0);
    weights.assign(q, 0.0);
    const int m = (q + 1) / 2;
    const double eps = 1e-15;
    for (int i = 0; i < m; ++i) {
        const double pi = std::acos(-1.0);
        double z = std::cos(pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(q) + 0.5));
        double z_prev = 0.0;
        double pp = 0.0;
        while (std::abs(z - z_prev) > eps) {
            double p1 = 1.0;
            double p2 = 0.0;
            for (int j = 1; j <= q; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = q * (z * p1 - p2) / (z * z - 1.0);
            z_prev = z;
            z = z_prev - p1 / pp;
        }
        const double x_left = -z;
        const double x_right = z;
        const double w = 2.0 / ((1.0 - z * z) * pp * pp);
        nodes[i] = 0.5 * (x_left + 1.0);
        nodes[q - 1 - i] = 0.5 * (x_right + 1.0);
        weights[i] = 0.5 * w;
        weights[q - 1 - i] = 0.5 * w;
    }
}

static inline double simplex_log_density_gauss(
    const std::vector<double>& h,
    const double volume,
    std::vector<double>& grad,
    const int q
) {
    const int n = h.size();
    const int d = n - 1;
    const double shift = *std::max_element(h.begin(), h.end());
    std::vector<double> nodes;
    std::vector<double> weights;
    gauss_legendre_unit_interval(q, nodes, weights);

    std::fill(grad.begin(), grad.end(), 0.0);
    double value_scaled = 0.0;
    std::vector<double> lambda(n, 0.0);

    if (d == 2) {
        for (int a = 0; a < q; ++a) {
            const double y1 = nodes[a];
            const double w1 = weights[a] * 2.0 * (1.0 - y1);
            for (int b = 0; b < q; ++b) {
                const double y2 = nodes[b];
                const double w = w1 * weights[b];
                lambda[0] = y1;
                lambda[1] = (1.0 - y1) * y2;
                lambda[2] = (1.0 - y1) * (1.0 - y2);
                double eta = 0.0;
                for (int i = 0; i < n; ++i) eta += lambda[i] * h[i];
                const double e = w * std::exp(eta - shift);
                value_scaled += e;
                for (int i = 0; i < n; ++i) grad[i] += e * lambda[i];
            }
        }
    } else if (d == 3) {
        for (int a = 0; a < q; ++a) {
            const double y1 = nodes[a];
            const double rem1 = 1.0 - y1;
            const double w1 = weights[a] * 3.0 * rem1 * rem1;
            for (int b = 0; b < q; ++b) {
                const double y2 = nodes[b];
                const double rem2 = 1.0 - y2;
                const double w2 = w1 * weights[b] * 2.0 * rem2;
                for (int c = 0; c < q; ++c) {
                    const double y3 = nodes[c];
                    const double w = w2 * weights[c];
                    lambda[0] = y1;
                    lambda[1] = rem1 * y2;
                    lambda[2] = rem1 * rem2 * y3;
                    lambda[3] = rem1 * rem2 * (1.0 - y3);
                    double eta = 0.0;
                    for (int i = 0; i < n; ++i) eta += lambda[i] * h[i];
                    const double e = w * std::exp(eta - shift);
                    value_scaled += e;
                    for (int i = 0; i < n; ++i) grad[i] += e * lambda[i];
                }
            }
        }
    } else {
        stop("Gauss fallback supports only triangles and tetrahedra.");
    }

    const double scale = volume * std::exp(shift);
    for (int i = 0; i < n; ++i) {
        grad[i] *= scale;
    }
    return scale * value_scaled;
}

static inline double simplex_log_density_integral_one(
    const std::vector<double>& h,
    const double volume,
    std::vector<double>& grad,
    const double taylor_radius,
    const double close_tol,
    const double taylor_tol,
    const int taylor_max_terms,
    const int gauss_order,
    int& method
) {
    const int n = h.size();
    const auto minmax = std::minmax_element(h.begin(), h.end());
    const double spread = *minmax.second - *minmax.first;
    if (spread <= taylor_radius) {
        method = 1;
        return simplex_log_density_taylor(h, volume, grad, taylor_tol, taylor_max_terms);
    }

    double min_gap = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            min_gap = std::min(min_gap, std::abs(h[i] - h[j]));
        }
    }
    if (min_gap <= close_tol * std::max(1.0, spread)) {
        method = 2;
        return simplex_log_density_gauss(h, volume, grad, gauss_order);
    }

    method = 0;
    return simplex_log_density_direct(h, volume, grad);
}

//' Integrate a log-linear density over simplexes
//'
//' For each row of `h`, computes `volume * E[exp(sum_i lambda_i h_i)]`
//' where `lambda` is uniform over the triangle or tetrahedron. The gradient
//' columns are the derivatives with respect to the corresponding vertex values.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List simplex_log_density_integral_cpp(
    const NumericMatrix& h,
    const NumericVector& volume,
    const double taylor_radius,
    const double close_tol,
    const double taylor_tol,
    const int taylor_max_terms,
    const int gauss_order,
    const int n_threads
) {
    const int n_simplex = h.nrow();
    const int n_vertex = h.ncol();
    if (n_vertex != 3 && n_vertex != 4) {
        stop("h must have three columns for triangles or four columns for tetrahedra.");
    }
    if (volume.size() != 1 && volume.size() != n_simplex) {
        stop("volume must have length 1 or one value per row of h.");
    }
    if (!R_finite(taylor_radius) || taylor_radius < 0.0) {
        stop("taylor_radius must be a non-negative finite scalar.");
    }
    if (!R_finite(close_tol) || close_tol < 0.0) {
        stop("close_tol must be a non-negative finite scalar.");
    }
    if (!R_finite(taylor_tol) || taylor_tol <= 0.0) {
        stop("taylor_tol must be a positive finite scalar.");
    }
    if (taylor_max_terms < 0) {
        stop("taylor_max_terms must be non-negative.");
    }
    if (gauss_order < 2) {
        stop("gauss_order must be at least 2.");
    }
    for (int i = 0; i < n_simplex; ++i) {
        for (int a = 0; a < n_vertex; ++a) {
            if (!R_finite(h(i, a))) {
                stop("h must contain only finite values.");
            }
        }
    }
    for (int i = 0; i < volume.size(); ++i) {
        if (!R_finite(volume[i]) || volume[i] < 0.0) {
            stop("volume values must be non-negative and finite.");
        }
    }

    NumericVector integral(n_simplex);
    NumericMatrix gradient(n_simplex, n_vertex);
    IntegerVector method(n_simplex);
    const int actual_threads = log_density_resolve_threads(n_threads);

#ifdef _OPENMP
#pragma omp parallel for num_threads(actual_threads)
#endif
    for (int i = 0; i < n_simplex; ++i) {
        std::vector<double> h_i(n_vertex);
        for (int a = 0; a < n_vertex; ++a) {
            h_i[a] = h(i, a);
        }
        const double vol_i = volume.size() == 1 ? volume[0] : volume[i];
        std::vector<double> grad_i(n_vertex, 0.0);
        int method_i = 0;
        const double value_i = simplex_log_density_integral_one(
            h_i,
            vol_i,
            grad_i,
            taylor_radius,
            close_tol,
            taylor_tol,
            taylor_max_terms,
            gauss_order,
            method_i
        );
        integral[i] = value_i;
        method[i] = method_i;
        for (int a = 0; a < n_vertex; ++a) {
            gradient(i, a) = grad_i[a];
        }
    }

    return List::create(
        _["integral"] = integral,
        _["gradient"] = gradient,
        _["method"] = method
    );
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
