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

struct BccPoint {
    int x;
    int y;
    int z;
};

static inline int sq_dist(const BccPoint& a, const BccPoint& b) {
    const int dx = a.x - b.x;
    const int dy = a.y - b.y;
    const int dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

static inline bool is_bcc_delaunay_tet(const BccPoint v[4]) {
    int d2[6];
    int idx = 0;

    for (int a = 0; a < 3; ++a) {
        for (int b = a + 1; b < 4; ++b) {
            d2[idx++] = sq_dist(v[a], v[b]);
        }
    }

    std::sort(d2, d2 + 6);

    return d2[0] == 3 && d2[1] == 3 && d2[2] == 3 &&
        d2[3] == 3 && d2[4] == 4 && d2[5] == 4;
}

static inline double det3(
    const double a11, const double a12, const double a13,
    const double a21, const double a22, const double a23,
    const double a31, const double a32, const double a33
) {
    return a11 * (a22 * a33 - a23 * a32) -
        a12 * (a21 * a33 - a23 * a31) +
        a13 * (a21 * a32 - a22 * a31);
}

static inline bool barycentric_weights(
    const BccPoint v[4],
    const double qx,
    const double qy,
    const double qz,
    double w[4],
    const double tol
) {
    const double a11 = static_cast<double>(v[0].x - v[3].x);
    const double a21 = static_cast<double>(v[0].y - v[3].y);
    const double a31 = static_cast<double>(v[0].z - v[3].z);

    const double a12 = static_cast<double>(v[1].x - v[3].x);
    const double a22 = static_cast<double>(v[1].y - v[3].y);
    const double a32 = static_cast<double>(v[1].z - v[3].z);

    const double a13 = static_cast<double>(v[2].x - v[3].x);
    const double a23 = static_cast<double>(v[2].y - v[3].y);
    const double a33 = static_cast<double>(v[2].z - v[3].z);

    const double b1 = qx - static_cast<double>(v[3].x);
    const double b2 = qy - static_cast<double>(v[3].y);
    const double b3 = qz - static_cast<double>(v[3].z);

    const double det_a = det3(
        a11, a12, a13,
        a21, a22, a23,
        a31, a32, a33
    );

    if (std::abs(det_a) <= tol) {
        return false;
    }

    w[0] = det3(
        b1, a12, a13,
        b2, a22, a23,
        b3, a32, a33
    ) / det_a;

    w[1] = det3(
        a11, b1, a13,
        a21, b2, a23,
        a31, b3, a33
    ) / det_a;

    w[2] = det3(
        a11, a12, b1,
        a21, a22, b2,
        a31, a32, b3
    ) / det_a;

    w[3] = 1.0 - w[0] - w[1] - w[2];

    for (int i = 0; i < 4; ++i) {
        if (w[i] < -tol || w[i] > 1.0 + tol) {
            return false;
        }
    }

    double sum_w = 0.0;
    for (int i = 0; i < 4; ++i) {
        if (std::abs(w[i]) <= tol) {
            w[i] = 0.0;
        } else if (std::abs(w[i] - 1.0) <= tol) {
            w[i] = 1.0;
        }
        sum_w += w[i];
    }

    if (sum_w != 0.0 && std::abs(sum_w - 1.0) > tol) {
        for (int i = 0; i < 4; ++i) {
            w[i] /= sum_w;
        }
    }

    return true;
}

static inline int lower_with_parity(const double x, const int parity) {
    return parity + 2 * static_cast<int>(std::floor((x - static_cast<double>(parity)) / 2.0));
}

// C++ backend for bcc_barycentric().
// [[Rcpp::export]]
List bcc_barycentric_cpp(
    const NumericMatrix& coords,
    const double s,
    const NumericVector& origin,
    const double tol = 1e-10,
    const int n_threads = 0
) {
    const int n = coords.nrow();

    if (coords.ncol() != 3) {
        stop("coords must have exactly 3 columns.");
    }
    if (!R_finite(s) || s <= 0.0) {
        stop("s must be a positive finite number.");
    }
    if (origin.size() != 3) {
        stop("origin must have length 3.");
    }
    for (int j = 0; j < 3; ++j) {
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
        if (
            !R_finite(coords(row, 0)) ||
            !R_finite(coords(row, 1)) ||
            !R_finite(coords(row, 2))
        ) {
            stop("coords contains non-finite values.");
        }
    }

    NumericMatrix weights(n, 4);
    NumericVector points(static_cast<R_xlen_t>(n) * 4 * 3);
    IntegerVector lattice(static_cast<R_xlen_t>(n) * 4 * 3);
    IntegerVector found_rows(n);

    points.attr("dim") = IntegerVector::create(n, 4, 3);
    lattice.attr("dim") = IntegerVector::create(n, 4, 3);

    const int actual_threads = resolve_threads(n_threads);

#ifdef _OPENMP
#pragma omp parallel for num_threads(actual_threads)
#endif
    for (int row = 0; row < n; ++row) {
        const double qx = (coords(row, 0) - origin[0]) / s;
        const double qy = (coords(row, 1) - origin[1]) / s;
        const double qz = (coords(row, 2) - origin[2]) / s;

        BccPoint candidates[16];
        int cidx = 0;

        for (int parity = 0; parity <= 1; ++parity) {
            const int x0 = lower_with_parity(qx, parity);
            const int y0 = lower_with_parity(qy, parity);
            const int z0 = lower_with_parity(qz, parity);

            const int xs[2] = {x0, x0 + 2};
            const int ys[2] = {y0, y0 + 2};
            const int zs[2] = {z0, z0 + 2};

            for (int ix = 0; ix < 2; ++ix) {
                for (int iy = 0; iy < 2; ++iy) {
                    for (int iz = 0; iz < 2; ++iz) {
                        candidates[cidx++] = BccPoint{xs[ix], ys[iy], zs[iz]};
                    }
                }
            }
        }

        bool found = false;
        BccPoint tet[4];
        double w[4] = {0.0, 0.0, 0.0, 0.0};

        for (int a = 0; a < 13 && !found; ++a) {
            tet[0] = candidates[a];
            for (int b = a + 1; b < 14 && !found; ++b) {
                tet[1] = candidates[b];
                for (int c = b + 1; c < 15 && !found; ++c) {
                    tet[2] = candidates[c];
                    for (int d = c + 1; d < 16 && !found; ++d) {
                        tet[3] = candidates[d];

                        if (!is_bcc_delaunay_tet(tet)) {
                            continue;
                        }

                        found = barycentric_weights(tet, qx, qy, qz, w, tol);
                    }
                }
            }
        }

        if (!found) {
            found_rows[row] = 0;
            continue;
        }
        found_rows[row] = 1;

        for (int a = 0; a < 4; ++a) {
            weights(row, a) = w[a];

            const BccPoint v = tet[a];
            const int normalized[3] = {v.x, v.y, v.z};

            for (int dim = 0; dim < 3; ++dim) {
                const R_xlen_t out_idx = row + static_cast<R_xlen_t>(n) * (a + 4 * dim);
                lattice[out_idx] = normalized[dim];
                points[out_idx] = origin[dim] + s * static_cast<double>(normalized[dim]);
            }
        }
    }

    for (int row = 0; row < n; ++row) {
        if (found_rows[row] == 0) {
            stop("No containing BCC tetrahedron found for row %d.", row + 1);
        }
    }

    colnames(weights) = CharacterVector::create("w1", "w2", "w3", "w4");

    return List::create(
        _["weights"] = weights,
        _["points"] = points,
        _["lattice"] = lattice
    );
}
