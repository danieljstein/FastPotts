#include <Rcpp.h>
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ------------------------------------------------------------
// nD Hilbert index from integer coordinates on [0, 2^bits - 1]
// dims must be 2 or 3 in this wrapper
// ------------------------------------------------------------
static inline uint64_t hilbert_index_nd(std::vector<uint32_t> x, int dims, int bits) {
    uint32_t M = 1u << (bits - 1);

    for (uint32_t Q = M; Q > 1; Q >>= 1) {
        uint32_t P = Q - 1;
        for (int i = 0; i < dims; i++) {
            if (x[i] & Q) {
                x[0] ^= P;
            } else {
                uint32_t t = (x[0] ^ x[i]) & P;
                x[0] ^= t;
                x[i] ^= t;
            }
        }
    }

    for (int i = 1; i < dims; i++) {
        x[i] ^= x[i - 1];
    }

    uint32_t t = 0;
    for (uint32_t Q = M; Q > 1; Q >>= 1) {
        if (x[dims - 1] & Q) {
            t ^= (Q - 1);
        }
    }
    for (int i = 0; i < dims; i++) {
        x[i] ^= t;
    }

    uint64_t h = 0;
    for (int b = bits - 1; b >= 0; --b) {
        for (int i = 0; i < dims; i++) {
            h = (h << 1) | ((x[i] >> b) & 1u);
        }
    }

    return h;
}

// ------------------------------------------------------------
// Scale raw coordinates to integer grid with shared aspect ratio
// ------------------------------------------------------------
static inline void scale_coords_global(
    const NumericMatrix& coords,
    int bits,
    IntegerMatrix& scaled,
    NumericVector& mins,
    double& global_extent
) {
    const int n = coords.nrow();
    const int dims = coords.ncol();
    const double grid_max = std::ldexp(1.0, bits) - 1.0; // 2^bits - 1

    mins = NumericVector(dims);

    for (int j = 0; j < dims; j++) {
        double mn = coords(0, j);
        if (!R_finite(mn)) stop("coords contains non-finite values.");
        for (int i = 1; i < n; i++) {
            double v = coords(i, j);
            if (!R_finite(v)) stop("coords contains non-finite values.");
            if (v < mn) mn = v;
        }
        mins[j] = mn;
    }

    NumericVector extents(dims);
    global_extent = 0.0;

    for (int j = 0; j < dims; j++) {
        double mx_shifted = 0.0;
        for (int i = 0; i < n; i++) {
            double shifted = coords(i, j) - mins[j];
            if (shifted > mx_shifted) mx_shifted = shifted;
        }
        extents[j] = mx_shifted;
        if (mx_shifted > global_extent) global_extent = mx_shifted;
    }

    if (global_extent == 0.0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < dims; j++) {
                scaled(i, j) = 0;
            }
        }
        return;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dims; j++) {
            double shifted = coords(i, j) - mins[j];
            double val = std::round((shifted / global_extent) * grid_max);

            if (val < 0.0) val = 0.0;
            if (val > grid_max) val = grid_max;

            scaled(i, j) = static_cast<int>(val);
        }
    }
}

// [[Rcpp::export]]
SEXP hilbert_index_scaled_cpp(NumericMatrix coords, int bits, bool return_scaled = false) {
    const int n = coords.nrow();
    const int dims = coords.ncol();

    if (dims != 2 && dims != 3) {
        stop("coords must have 2 or 3 columns.");
    }
    if (n < 1) {
        stop("coords must have at least one row.");
    }
    if (bits < 1 || bits > 21) {
        stop("bits must be between 1 and 21 for this integer-scaled wrapper.");
    }

    IntegerMatrix scaled(n, dims);
    NumericVector mins(dims);
    double global_extent = 0.0;

    scale_coords_global(coords, bits, scaled, mins, global_extent);

    NumericVector index(n);
    std::vector<uint32_t> x(dims);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dims; j++) {
            x[j] = static_cast<uint32_t>(scaled(i, j));
        }
        uint64_t h = hilbert_index_nd(x, dims, bits);
        index[i] = static_cast<double>(h);
    }

    if (!return_scaled) {
        return index;
    }

    return List::create(
        _["index"] = index,
        _["scaled_coords"] = scaled,
        _["mins"] = mins,
        _["global_extent"] = global_extent,
        _["bits"] = bits
    );
}
