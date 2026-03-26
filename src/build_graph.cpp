// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <utility>
using namespace Rcpp;

/*
 * Binary-search for a row index inside one CSC column.
 *
 * In a dgCMatrix, column `col` occupies indices p[col]:(p[col+1]-1),
 * and row indices i[...] are sorted within each column.
 *
 * Returns the nonzero index if found, or -1 otherwise.
 */
static inline int find_in_col(
    const IntegerVector& p,
    const IntegerVector& i,
    const int col,
    const int target_row
) {
    int left = p[col];
    int right = p[col + 1] - 1;

    while (left <= right) {
        const int mid = left + (right - left) / 2;
        const int val = i[mid];

        if (val == target_row) return mid;
        if (val < target_row) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return -1;
}

static inline std::uint64_t edge_key(const int from, const int to) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(from)) << 32) |
        static_cast<std::uint32_t>(to);
}

//' Build Potts-LBP graph inputs from a sparse adjacency matrix
//'
//' Converts a square sparse adjacency matrix of class `dgCMatrix` into the
//' CSR-like graph representation expected by `potts_lbp()`.
//'
//' The input matrix is interpreted as a directed adjacency matrix stored in
//' CSC format. If the graph is undirected, each edge should appear twice:
//' once at `(i, j)` and once at `(j, i)`. The function constructs:
//'
//' - `adj_ptr`: CSR row pointer array
//' - `adj_idx`: destination node index for each directed edge
//' - `rev_idx`: index of the reverse directed edge
//' - `edge_weights`: edge weights in the same order as `adj_idx`
//'
//' Node indices in the returned object are **0-based**, matching the C++
//' indexing expected by `potts_lbp()`.
//'
//' @param mat A square sparse matrix of class `dgCMatrix`. Off-diagonal nonzero
//'   entries define directed edges; diagonal entries are ignored by default.
//' @param drop_diagonal Logical; if `TRUE` (default), diagonal entries are
//'   ignored. If `FALSE`, diagonal entries are retained, but then each diagonal
//'   entry is its own reverse.
//' @param check_reverse Logical; if `TRUE` (default), each directed edge must
//'   have a matching reverse edge or an error is thrown.
//'
//' @return A list with components:
//' \describe{
//'   \item{adj_ptr}{Integer vector of length `N + 1` giving CSR row pointers.}
//'   \item{adj_idx}{Integer vector of length `E` of 0-based destination nodes.}
//'   \item{rev_idx}{Integer vector of length `E`; `rev_idx[e]` is the index of
//'     the reverse directed edge for edge `e`.}
//'   \item{edge_weights}{Numeric vector of length `E` of edge weights in the
//'     same order as `adj_idx`.}
//'   \item{n_nodes}{Number of nodes.}
//'   \item{n_edges}{Number of directed edges after optional diagonal removal.}
//' }
//'
//' @details
//' A `dgCMatrix` is stored in compressed sparse column (CSC) format. This helper
//' first counts outgoing edges per row in order to build a CSR representation,
//' then fills the directed edge list, and finally computes reverse-edge indices.
//'
//' Reverse edges are found efficiently by binary-searching for entry `(j, i)`
//' in column `i` for each entry `(i, j)`. For sparse k-NN graphs, this is
//' typically very fast because each column contains only a small number of
//' nonzeros.
//'
//' The returned edge order is by source row, then by the original CSC scan order.
//' This is suitable for the LBP implementation, which only requires consistency
//' among `adj_ptr`, `adj_idx`, `rev_idx`, and `edge_weights`.
//'
//' @examples
//' \dontrun{
//' library(Matrix)
//'
//' W <- sparseMatrix(
//'   i = c(1,2,1,3,2,3),
//'   j = c(2,1,3,1,3,2),
//'   x = c(0.7,0.7,0.4,0.4,0.9,0.9),
//'   dims = c(3,3)
//' )
//'
//' g <- build_potts_lbp_graph(W)
//'
//' str(g)
//' # g$adj_ptr
//' # g$adj_idx
//' # g$rev_idx
//' # g$edge_weights
//' }
//'
//' @export
// [[Rcpp::export]]
List build_potts_lbp_graph(
    SEXP mat,
    const bool drop_diagonal = true,
    const bool check_reverse = true
) {
    S4 A(mat);

    if (!A.is("dgCMatrix")) {
        stop("`mat` must be a dgCMatrix.");
    }

    IntegerVector p = A.slot("p");
    IntegerVector i = A.slot("i");
    NumericVector x = A.slot("x");
    IntegerVector Dim = A.slot("Dim");

    const int N = Dim[0];
    const int M = Dim[1];

    if (N != M) {
        stop("`mat` must be square.");
    }

    const int nnz = i.size();

    // Step 1: count outgoing edges per row after optional diagonal removal
    IntegerVector row_count(N);
    for (int col = 0; col < N; ++col) {
        for (int idx = p[col]; idx < p[col + 1]; ++idx) {
            const int row = i[idx];
            if (drop_diagonal && row == col) continue;
            row_count[row] += 1;
        }
    }

    // Step 2: build CSR row pointers
    IntegerVector adj_ptr(N + 1);
    adj_ptr[0] = 0;
    for (int r = 0; r < N; ++r) {
        adj_ptr[r + 1] = adj_ptr[r] + row_count[r];
    }

    const int E = adj_ptr[N];

    IntegerVector adj_idx(E);
    NumericVector edge_weights(E);

    // Temporary write cursor per row
    IntegerVector next_pos(N);
    for (int r = 0; r < N; ++r) {
        next_pos[r] = adj_ptr[r];
    }

    // We also keep track of which original CSC entry generated each CSR edge.
    // This lets us compute reverse-edge indices efficiently afterward.
    IntegerVector csr_to_csc(E);
    IntegerVector csc_to_csr(nnz, NA_INTEGER);

    // Step 3: fill CSR arrays by scanning CSC entries
    for (int col = 0; col < N; ++col) {
        for (int idx = p[col]; idx < p[col + 1]; ++idx) {
            const int row = i[idx];
            if (drop_diagonal && row == col) continue;

            const int pos = next_pos[row]++;
            adj_idx[pos] = col;          // 0-based destination node
            edge_weights[pos] = x[idx];
            csr_to_csc[pos] = idx;
            csc_to_csr[idx] = pos;
        }
    }

    // Step 4: compute reverse-edge indices
    IntegerVector rev_idx(E);

    for (int edge = 0; edge < E; ++edge) {
        const int idx = csr_to_csc[edge];
        int col = -1;

        // Recover the CSC column containing idx.
        // Because p is monotone, binary search is appropriate.
        {
            int left = 0;
            int right = N - 1;
            while (left <= right) {
                const int mid = left + (right - left) / 2;
                if (idx < p[mid]) {
                    right = mid - 1;
                } else if (idx >= p[mid + 1]) {
                    left = mid + 1;
                } else {
                    col = mid;
                    break;
                }
            }
        }

        if (col < 0) {
            stop("Internal error: could not recover CSC column for edge %d.", edge + 1);
        }

        const int row = i[idx];

        // Current edge corresponds to (row -> col) in CSR terms,
        // i.e. matrix entry at (row, col). Its reverse is (col, row),
        // which must be found in CSC column = row, with row index = col.
        int ridx;
        if (!drop_diagonal && row == col) {
            ridx = idx;  // diagonal entry is its own reverse
        } else {
            ridx = find_in_col(p, i, row, col);
        }

        if (ridx < 0) {
            if (check_reverse) {
                stop(
                    "Reverse edge not found for matrix entry (%d, %d).",
                    row + 1, col + 1
                );
            } else {
                rev_idx[edge] = NA_INTEGER;
                continue;
            }
        }

        const int rev_edge = csc_to_csr[ridx];
        if (rev_edge == NA_INTEGER) {
            if (check_reverse) {
                stop(
                    "Reverse edge exists structurally but was dropped for entry (%d, %d).",
                    row + 1, col + 1
                );
            } else {
                rev_idx[edge] = NA_INTEGER;
                continue;
            }
        }

        rev_idx[edge] = rev_edge;
    }

    return List::create(
        _["adj_ptr"] = adj_ptr,
        _["adj_idx"] = adj_idx,
        _["rev_idx"] = rev_idx,
        _["edge_weights"] = edge_weights,
        _["n_nodes"] = N,
        _["n_edges"] = E
    );
}

//' Build Potts-LBP graph inputs directly from an edge list
//'
//' Builds the CSR-like graph representation expected by `potts_lbp()` from a
//' directed edge list without first materializing a sparse matrix.
//'
//' Duplicate directed edges are deduplicated so that each directed edge appears
//' at most once in the returned graph. If duplicate `(from, to)` entries are
//' supplied with different weights, an error is thrown.
//'
//' Node indices supplied in `from` and `to` are **1-based R indices**. Returned
//' graph indices remain **0-based**, matching the C++ indexing expected by
//' `potts_lbp()`.
//'
//' @param from Integer vector of 1-based source node indices.
//' @param to Integer vector of 1-based destination node indices.
//' @param weights Numeric vector of edge weights. Must either have length 1,
//'   in which case the same weight is used for every input edge, or the same
//'   length as `from`.
//' @param n_nodes Integer number of nodes in the graph.
//' @param symmetric Logical; if `TRUE`, each input edge `(u, v)` is treated as
//'   undirected and inserted in both directions `(u, v)` and `(v, u)`.
//' @param check_reverse Logical; if `TRUE` (default), each directed edge must
//'   have a matching reverse edge after duplicate aggregation.
//'
//' @return A list with the same structure as [build_potts_lbp_graph()].
//'
//' @keywords internal
// [[Rcpp::export]]
List build_potts_lbp_graph_from_edges(
    const IntegerVector& from,
    const IntegerVector& to,
    const NumericVector& weights,
    const int n_nodes,
    const bool symmetric = false,
    const bool check_reverse = true
) {
    const R_xlen_t n_input_edges = from.size();

    if (to.size() != n_input_edges) {
        stop("`from` and `to` must have the same length.");
    }
    if (n_nodes <= 0) {
        stop("`n_nodes` must be positive.");
    }
    if (!(weights.size() == 1 || weights.size() == n_input_edges)) {
        stop("`weights` must have length 1 or the same length as `from`.");
    }

    std::unordered_map<std::uint64_t, double> edge_weights_map;
    edge_weights_map.reserve(static_cast<std::size_t>(symmetric ? 2 * n_input_edges : n_input_edges));

    auto insert_edge = [&edge_weights_map](const int src0, const int dst0, const double weight) {
        const std::uint64_t key = edge_key(src0, dst0);
        const auto existing = edge_weights_map.find(key);
        if (existing == edge_weights_map.end()) {
            edge_weights_map.emplace(key, weight);
        } else if (existing->second != weight) {
            stop(
                "Duplicate edge (%d, %d) supplied with conflicting weights.",
                src0 + 1, dst0 + 1
            );
        }
    };

    for (R_xlen_t idx = 0; idx < n_input_edges; ++idx) {
        const int src = from[idx];
        const int dst = to[idx];

        if (src == NA_INTEGER || dst == NA_INTEGER) {
            stop("`from` and `to` must not contain missing values.");
        }
        if (src < 1 || src > n_nodes || dst < 1 || dst > n_nodes) {
            stop("Edge endpoints must be between 1 and n_nodes.");
        }

        const double weight = weights.size() == 1 ? weights[0] : weights[idx];
        if (!R_finite(weight)) {
            stop("`weights` must be finite.");
        }

        insert_edge(src - 1, dst - 1, weight);
        if (symmetric && src != dst) {
            insert_edge(dst - 1, src - 1, weight);
        }
    }

    std::vector<int> row_count(n_nodes, 0);
    for (const auto& kv : edge_weights_map) {
        const int src = static_cast<int>(kv.first >> 32);
        row_count[src] += 1;
    }

    IntegerVector adj_ptr(n_nodes + 1);
    adj_ptr[0] = 0;
    for (int node = 0; node < n_nodes; ++node) {
        adj_ptr[node + 1] = adj_ptr[node] + row_count[node];
    }

    const int n_edges = adj_ptr[n_nodes];
    IntegerVector adj_idx(n_edges);
    NumericVector edge_weights_out(n_edges);
    IntegerVector rev_idx(n_edges);

    std::vector<std::vector<std::pair<int, double>>> rows(n_nodes);
    for (int node = 0; node < n_nodes; ++node) {
        rows[node].reserve(row_count[node]);
    }

    for (const auto& kv : edge_weights_map) {
        const int src = static_cast<int>(kv.first >> 32);
        const int dst = static_cast<int>(kv.first & 0xffffffffu);
        rows[src].push_back(std::make_pair(dst, kv.second));
    }

    std::unordered_map<std::uint64_t, int> edge_position;
    edge_position.reserve(static_cast<std::size_t>(n_edges));

    int edge = 0;
    for (int src = 0; src < n_nodes; ++src) {
        auto& row = rows[src];
        std::sort(
            row.begin(),
            row.end(),
            [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs) {
                return lhs.first < rhs.first;
            }
        );

        for (const auto& entry : row) {
            adj_idx[edge] = entry.first;
            edge_weights_out[edge] = entry.second;
            edge_position[edge_key(src, entry.first)] = edge;
            edge += 1;
        }
    }

    for (int src = 0; src < n_nodes; ++src) {
        for (int pos = adj_ptr[src]; pos < adj_ptr[src + 1]; ++pos) {
            const int dst = adj_idx[pos];
            const auto rev = edge_position.find(edge_key(dst, src));
            if (rev == edge_position.end()) {
                if (check_reverse) {
                    stop(
                        "Reverse edge not found for edge (%d, %d).",
                        src + 1, dst + 1
                    );
                }
                rev_idx[pos] = NA_INTEGER;
            } else {
                rev_idx[pos] = rev->second;
            }
        }
    }

    return List::create(
        _["adj_ptr"] = adj_ptr,
        _["adj_idx"] = adj_idx,
        _["rev_idx"] = rev_idx,
        _["edge_weights"] = edge_weights_out,
        _["n_nodes"] = n_nodes,
        _["n_edges"] = n_edges
    );
}
