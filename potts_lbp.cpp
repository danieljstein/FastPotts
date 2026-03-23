// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

/*
 * Numerically stable log-sum-exp for a length-K float vector.
 *
 * Given x[0:(K-1)], returns
 *
 *   log( sum_k exp(x[k]) )
 *
 * using the standard "subtract the maximum" stabilization trick.
 *
 * Parameters
 * ----------
 * x : pointer to the first element of a length-K float array
 * K : number of entries
 *
 * Returns
 * -------
 * A float containing log(sum(exp(x))).
 */
static inline float logsumexp_float(const float* x, int K) {
    float mx = x[0];
    for (int i = 1; i < K; ++i) {
        if (x[i] > mx) mx = x[i];
    }

    double s = 0.0;
    for (int i = 0; i < K; ++i) {
        s += std::exp((double)x[i] - (double)mx);
    }
    return mx + (float)std::log(s);
}

//' Loopy belief propagation for a Potts model on a sparse directed graph
//'
//' Runs sum-product loopy belief propagation (LBP) for a Potts model with
//' node-specific unary log-potentials and edge-specific attractive pairwise
//' log-potentials.
//'
//' The graph is supplied in CSR-like directed adjacency form. Each undirected
//' edge should appear twice in `adj_idx`: once as `u -> v` and once as
//' `v -> u`. The vector `rev_idx` maps each directed edge to the index of its
//' reverse directed edge.
//'
//' The pairwise potential is parameterized as
//' \deqn{
//' \psi_{uv}(a,b) = \exp(w_{uv}) \text{ if } a=b,\quad 1 \text{ otherwise},
//' }
//' where `edge_weights[edge] = w_{uv}` is the log same-label preference for the
//' directed edge `u -> v`. In most applications the two directions of an
//' undirected edge will have the same weight.
//'
//' The unary potential is parameterized by `node_potential`, whose `(i, k)` entry
//' is the log unary potential for node `i` and label `k`.
//'
//' Message updates are performed in log-space for numerical stability. The
//' returned node marginals are normalized probabilities on the original scale.
//'
//' @param adj_ptr Integer vector of length `N + 1` giving CSR row pointers for
//'   the directed adjacency structure. Directed neighbors of node `u` are stored
//'   in `adj_idx[adj_ptr[u] : adj_ptr[u + 1] - 1]` using 0-based C++ indexing.
//' @param adj_idx Integer vector of length `E` containing the destination node
//'   for each directed edge. Values are expected to be 0-based node indices.
//' @param rev_idx Integer vector of length `E` such that `rev_idx[e]` is the
//'   directed-edge index of the reverse of edge `e`.
//' @param edge_weights Numeric vector of length `E` containing log pairwise
//'   Potts weights for each directed edge.
//' @param node_potential Numeric matrix of dimension `N x K` containing log
//'   unary potentials.
//' @param max_iter Maximum number of LBP sweeps.
//' @param damping Damping parameter in `(0, 1]`. Values below 1 blend each new
//'   message with the previous message to improve stability on loopy graphs.
//' @param tol Convergence tolerance. Iteration stops early when the maximum
//'   absolute message change in a sweep is below `tol`.
//' @param synchronous Logical; if `TRUE`, uses synchronous updates with a second
//'   message buffer. If `FALSE` (default), uses in-place asynchronous updates.
//'
//' @return A list with components:
//' \describe{
//'   \item{marginals}{`N x K` numeric matrix of approximate node marginals.}
//'   \item{iterations}{Number of sweeps performed.}
//'   \item{max_delta}{Maximum absolute message change in the final sweep.}
//' }
//'
//' @details
//' Let
//' \deqn{
//' H_i(k) = \log \phi_i(k) + \sum_{\ell \in N(i)} \log m_{\ell \to i}(k).
//' }
//' For the directed edge `u -> v`, the cavity field is
//' \deqn{
//' h_{u \to v}(k) = H_u(k) - \log m_{v \to u}(k).
//' }
//' With the Potts parameterization above, the message update simplifies to
//' \deqn{
//' m_{u \to v}(k) \propto S_{u \to v} +
//'   (\exp(w_{uv}) - 1)\exp(h_{u \to v}(k)),
//' }
//' where
//' \deqn{
//' S_{u \to v} = \sum_a \exp(h_{u \to v}(a)).
//' }
//' This reduces each per-edge update from `O(K^2)` for a general pairwise model
//' to `O(K)` for the Potts model.
//'
//' @examples
//' \dontrun{
//' # Suppose graph contains directed edges for both directions of each
//' # undirected edge, stored in CSR form.
//' fit <- potts_lbp(
//'   adj_ptr = adj_ptr,
//'   adj_idx = adj_idx,
//'   rev_idx = rev_idx,
//'   edge_weights = edge_weights,
//'   node_potential = node_potential,
//'   max_iter = 50,
//'   damping = 0.5,
//'   tol = 1e-4,
//'   synchronous = FALSE
//' )
//'
//' fit$marginals
//' fit$iterations
//' fit$max_delta
//' }
//'
//' @export
// [[Rcpp::export]]
List potts_lbp(
    const IntegerVector& adj_ptr,
    const IntegerVector& adj_idx,
    const IntegerVector& rev_idx,
    const NumericVector& edge_weights,     // log of pairwise potentials
    const NumericMatrix& node_potential,   // log of unary potentials
    const int max_iter = 50,
    const float damping = 0.5,
    const float tol = 1e-4,
    const bool synchronous = false
) {
    const int N = node_potential.nrow();   // number of nodes
    const int K = node_potential.ncol();   // number of labels
    const int E = adj_idx.size();          // number of directed edges

    if (adj_ptr.size() != N + 1) stop("adj_ptr must have length N+1");
    if (rev_idx.size() != E || edge_weights.size() != E) stop("edge vectors must match");
    if (damping <= 0.0 || damping > 1.0) stop("damping must be in (0,1]");

    const float one_minus_damping = 1.0f - damping;

    // Edge weights stored internally as float:
    // w[e] is the log same-label preference for directed edge e.
    std::vector<float> w(E);
    for (int e = 0; e < E; ++e) {
        w[e] = (float) edge_weights[e];
    }

    // Log-messages for directed edges.
    // Message for directed edge e and label k is stored at msg[e * K + k].
    std::vector<float> msg((size_t)E * K, 0.0f);
    std::vector<float> msg_new;
    if (synchronous) {
        msg_new.assign((size_t)E * K, 0.0f);
    }

    // H_i(k) = log phi_i(k) + sum_{l in N(i)} log m_{l->i}(k)
    // This caches the full incoming field at each node/label combination.
    std::vector<float> H((size_t)N * K);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < K; ++k) {
            H[(size_t)i * K + k] = (float)node_potential(i, k);
        }
    }

    // Temporary storage for one edge update.
    std::vector<float> tmp_h(K);
    std::vector<float> tmp_msg(K);

    int iters = 0;
    float max_delta = 0.0f;

    for (int iter = 0; iter < max_iter; ++iter) {
        max_delta = 0.0f;

        for (int u = 0; u < N; ++u) {
            for (int edge = adj_ptr[u]; edge < adj_ptr[u + 1]; ++edge) {
                const int v = adj_idx[edge];
                const int rev = rev_idx[edge];

                const size_t u_offset = (size_t)u * K;
                const size_t v_offset = (size_t)v * K;
                const size_t edge_offset = (size_t)edge * K;
                const size_t rev_offset = (size_t)rev * K;

                // Cavity field for message u -> v:
                // h_{u->v}(k) = H_u(k) - log m_{v->u}(k)
                for (int k = 0; k < K; ++k) {
                    tmp_h[k] = H[u_offset + k] - msg[rev_offset + k];
                }

                // log S where S = sum_a exp(h[a])
                const float logS = logsumexp_float(tmp_h.data(), K);

                // alpha = exp(w_uv) - 1 for the simplified Potts potential
                const double alpha = std::exp((double)w[edge]) - 1.0;

                // log m(k) = log( S + (exp(w)-1) * exp(h[k]) )
                for (int k = 0; k < K; ++k) {
                    double val;
                    if (alpha <= 0.0) {
                        val = (double)logS;
                    } else {
                        const double a0 = (double)logS;
                        const double a1 = std::log(alpha) + (double)tmp_h[k];
                        const double mx = (a0 > a1 ? a0 : a1);
                        val = mx + std::log(std::exp(a0 - mx) + std::exp(a1 - mx));
                    }
                    tmp_msg[k] = (float)val;
                }

                // Normalize message in log-space so that exp(message) sums to 1.
                const float lse = logsumexp_float(tmp_msg.data(), K);
                for (int k = 0; k < K; ++k) {
                    tmp_msg[k] -= lse;
                }

                if (synchronous) {
                    // Compute next-iteration messages into msg_new.
                    for (int k = 0; k < K; ++k) {
                        const float oldv = msg[edge_offset + k];
                        const float newv = one_minus_damping * oldv + damping * tmp_msg[k];
                        msg_new[edge_offset + k] = newv;

                        const float d = std::fabs(newv - oldv);
                        if (d > max_delta) max_delta = d;
                    }
                } else {
                    // In-place asynchronous update:
                    // update message u -> v and immediately fold its change into H[v].
                    for (int k = 0; k < K; ++k) {
                        const float oldv = msg[edge_offset + k];
                        const float newv = one_minus_damping * oldv + damping * tmp_msg[k];
                        msg[edge_offset + k] = newv;

                        // Message u -> v contributes to the incoming field of v.
                        H[v_offset + k] += (newv - oldv);

                        const float d = std::fabs(newv - oldv);
                        if (d > max_delta) max_delta = d;
                    }
                }
            }
        }

        if (synchronous) {
            // Rebuild H from scratch from the new message buffer.
            for (int i = 0; i < N; ++i) {
                for (int k = 0; k < K; ++k) {
                    H[(size_t)i * K + k] = (float)node_potential(i, k);
                }
            }

            for (int u = 0; u < N; ++u) {
                for (int edge = adj_ptr[u]; edge < adj_ptr[u + 1]; ++edge) {
                    const int v = adj_idx[edge];
                    const size_t edge_offset = (size_t)edge * K;
                    const size_t v_offset = (size_t)v * K;
                    for (int k = 0; k < K; ++k) {
                        H[v_offset + k] += msg_new[edge_offset + k];
                    }
                }
            }

            msg.swap(msg_new);
        }

        iters = iter + 1;
        if (max_delta < tol) break;
    }

    // Convert final log-beliefs H into normalized node marginals.
    NumericMatrix marginals(N, K);
    for (int i = 0; i < N; ++i) {
        const size_t offset = (size_t)i * K;
        const float lse = logsumexp_float(&H[offset], K);
        for (int k = 0; k < K; ++k) {
            marginals(i, k) = std::exp((double)(H[offset + k] - lse));
        }
    }

    return List::create(
        _["marginals"] = marginals,
        _["iterations"] = iters,
        _["max_delta"] = (double)max_delta
    );
}
