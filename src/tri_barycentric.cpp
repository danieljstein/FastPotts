// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <cmath>
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

// C++ backend for tri_barycentric().
// [[Rcpp::export]]
List tri_barycentric_cpp(
    const NumericMatrix& coords,
    const double s,
    const NumericVector& origin,
    const double tol = 1e-10,
    const int n_threads = 0
) {
    const int n = coords.nrow();
    const double sqrt3 = std::sqrt(3.0);

    if (coords.ncol() != 2) {
        stop("coords must have exactly 2 columns.");
    }
    if (!R_finite(s) || s <= 0.0) {
        stop("s must be a positive finite number.");
    }
    if (origin.size() != 2) {
        stop("origin must have length 2.");
    }
    for (int j = 0; j < 2; ++j) {
        if (!R_finite(origin[j])) {
            stop("origin must contain only finite values.");
        }
    }
    if (!R_finite(tol) || tol < 0.0) {
        stop("tol must be a non-negative finite number.");
    }
    if (n_threads < 0) {
        stop("n_threads must be NULL or a positive integer.");
    }
    for (int row = 0; row < n; ++row) {
        if (!R_finite(coords(row, 0)) || !R_finite(coords(row, 1))) {
            stop("coords contains non-finite values.");
        }
    }

    NumericMatrix weights(n, 3);
    NumericVector points(static_cast<R_xlen_t>(n) * 3 * 2);
    IntegerVector lattice(static_cast<R_xlen_t>(n) * 3 * 2);

    points.attr("dim") = IntegerVector::create(n, 3, 2);
    lattice.attr("dim") = IntegerVector::create(n, 3, 2);

    const int actual_threads = resolve_threads(n_threads);

#ifdef _OPENMP
#pragma omp parallel for num_threads(actual_threads)
#endif
    for (int row = 0; row < n; ++row) {
        const double x = coords(row, 0) - origin[0];
        const double y = coords(row, 1) - origin[1];

        const double v = 2.0 * y / (sqrt3 * s);
        const double u = x / s - 0.5 * v;

        const int i = static_cast<int>(std::floor(u));
        const int j = static_cast<int>(std::floor(v));
        double r = u - static_cast<double>(i);
        double t = v - static_cast<double>(j);

        if (std::abs(r) <= tol) r = 0.0;
        if (std::abs(t) <= tol) t = 0.0;
        if (std::abs(r - 1.0) <= tol) r = 1.0;
        if (std::abs(t - 1.0) <= tol) t = 1.0;

        int ij[3][2];
        double w[3];

        if (r + t <= 1.0 + tol) {
            ij[0][0] = i;
            ij[0][1] = j;
            ij[1][0] = i + 1;
            ij[1][1] = j;
            ij[2][0] = i;
            ij[2][1] = j + 1;

            w[0] = 1.0 - r - t;
            w[1] = r;
            w[2] = t;
        } else {
            ij[0][0] = i + 1;
            ij[0][1] = j + 1;
            ij[1][0] = i;
            ij[1][1] = j + 1;
            ij[2][0] = i + 1;
            ij[2][1] = j;

            w[0] = r + t - 1.0;
            w[1] = 1.0 - r;
            w[2] = 1.0 - t;
        }

        double sum_w = 0.0;
        for (int a = 0; a < 3; ++a) {
            if (std::abs(w[a]) <= tol) {
                w[a] = 0.0;
            } else if (std::abs(w[a] - 1.0) <= tol) {
                w[a] = 1.0;
            }
            sum_w += w[a];
        }
        if (sum_w != 0.0 && std::abs(sum_w - 1.0) > tol) {
            for (int a = 0; a < 3; ++a) {
                w[a] /= sum_w;
            }
        }

        for (int a = 0; a < 3; ++a) {
            weights(row, a) = w[a];

            const int li = ij[a][0];
            const int lj = ij[a][1];
            const double px = origin[0] + s * (static_cast<double>(li) + 0.5 * static_cast<double>(lj));
            const double py = origin[1] + s * (sqrt3 / 2.0) * static_cast<double>(lj);

            const R_xlen_t x_idx = row + static_cast<R_xlen_t>(n) * (a + 3 * 0);
            const R_xlen_t y_idx = row + static_cast<R_xlen_t>(n) * (a + 3 * 1);

            lattice[x_idx] = li;
            lattice[y_idx] = lj;
            points[x_idx] = px;
            points[y_idx] = py;
        }
    }

    colnames(weights) = CharacterVector::create("w1", "w2", "w3");

    return List::create(
        _["weights"] = weights,
        _["points"] = points,
        _["lattice"] = lattice
    );
}
