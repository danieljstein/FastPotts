// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

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

    NumericVector grad(par.size());
    std::vector<double> f(n_cell_types);
    std::vector<double> log_post(n_cell_types);
    std::vector<double> p(n_cell_types);
    std::vector<double> q(n_cell_types);

    double objective = 0.0;

    for (int i = 0; i < n; ++i) {
        std::fill(f.begin(), f.end(), 0.0);

        for (int a = 0; a < n_active; ++a) {
            const int m = basis_id(i, a);
            const double phi = basis_weight(i, a);

            if (m < 0 || m >= n_basis) {
                stop("basis_id contains an out-of-range basis index.");
            }

            for (int k = 0; k < n_free; ++k) {
                f[k] += phi * par[m + n_basis * k];
            }
        }

        const double log_z_prior = log_sum_exp(f);
        for (int k = 0; k < n_cell_types; ++k) {
            p[k] = std::exp(f[k] - log_z_prior);
        }

        const int g = gene_index[i];
        if (g < 0 || g >= log_signature.nrow()) {
            stop("gene_index contains an out-of-range gene index.");
        }

        for (int k = 0; k < n_cell_types; ++k) {
            log_post[k] = log_signature(g, k) + f[k];
        }

        const double log_z_post = log_sum_exp(log_post);
        objective -= log_z_post - log_z_prior;

        for (int k = 0; k < n_cell_types; ++k) {
            q[k] = std::exp(log_post[k] - log_z_post);
        }

        for (int a = 0; a < n_active; ++a) {
            const int m = basis_id(i, a);
            const double phi = basis_weight(i, a);

            for (int k = 0; k < n_free; ++k) {
                grad[m + n_basis * k] += phi * (p[k] - q[k]);
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

                objective += 0.5 * lambda * diff * diff;
                grad[idx1] += lambda * diff;
                grad[idx2] -= lambda * diff;
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
    const int n_cell_types
) {
    const int n = basis_id.nrow();
    const int n_active = basis_id.ncol();
    const int n_free = n_cell_types - 1;

    if (par.size() != n_basis * n_free) {
        stop("par has incompatible length.");
    }

    NumericMatrix prior(n, n_cell_types);
    NumericMatrix posterior(n, n_cell_types);
    NumericMatrix logits(n, n_cell_types);

    std::vector<double> f(n_cell_types);
    std::vector<double> log_post(n_cell_types);

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

    return List::create(
        _["prior"] = prior,
        _["posterior"] = posterior,
        _["logits"] = logits
    );
}
