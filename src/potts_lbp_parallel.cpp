// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

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

//' Parallel synchronous loopy belief propagation for a Potts model
//'
//' Synchronous sum-product loopy belief propagation (LBP) for a Potts model on
//' a sparse directed graph stored in CSR-like form.
//'
//' This version uses OpenMP to parallelize the message-update sweep over source
//' nodes. The rebuild of the cached incoming fields is performed serially.
//'
//' @param adj_ptr Integer vector of length `N + 1` giving CSR row pointers.
//' @param adj_idx Integer vector of length `E` containing 0-based destination
//'   node indices for each directed edge.
//' @param rev_idx Integer vector of length `E` such that `rev_idx[e]` is the
//'   directed-edge index of the reverse of edge `e`.
//' @param edge_weights Numeric vector of length `E` containing log same-label
//'   Potts weights for each directed edge.
//' @param node_potential Numeric matrix of dimension `N x K` containing log
//'   unary potentials.
//' @param max_iter Maximum number of LBP sweeps.
//' @param damping Damping parameter in `(0, 1]`. Default is `1.0`.
//' @param tol Convergence tolerance on the maximum absolute message change per
//'   sweep. Default is `1e-2`.
//' @param n_threads Number of OpenMP threads to use. If `NULL`, uses the
//'   current OpenMP default.
//' @param verbose If `TRUE`, prints the maximum and mean message change at each iteration.
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
//' @export
// [[Rcpp::export]]
List potts_lbp_parallel_cpp(
    const IntegerVector& adj_ptr,
    const IntegerVector& adj_idx,
    const IntegerVector& rev_idx,
    const NumericVector& edge_weights,     // log same-label Potts weights
    const NumericMatrix& node_potential,   // log unary potentials
    const int max_iter = 50,
    const float damping = 1.0,
    const float tol = 1e-2,
    Nullable<int> n_threads = R_NilValue,
    const bool verbose = false
) {
    const int N = node_potential.nrow();   // number of nodes
    const int K = node_potential.ncol();   // number of labels
    const int E = adj_idx.size();          // number of directed edges

    if (adj_ptr.size() != N + 1) stop("adj_ptr must have length N+1");
    if (rev_idx.size() != E || edge_weights.size() != E) stop("edge vectors must match");
    if (damping <= 0.0f || damping > 1.0f) stop("damping must be in (0,1]");
    if (n_threads.isNotNull()) {
        int nt = as<int>(n_threads);
        if (nt <= 0) stop("n_threads must be positive");
    }

#ifdef _OPENMP
    if (n_threads.isNotNull()) {
        omp_set_num_threads(as<int>(n_threads));
    }
#endif

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
    std::vector<float> msg_new((size_t)E * K, 0.0f);

    // H_i(k) = log phi_i(k) + sum_{l in N(i)} log m_{l->i}(k)
    // This caches the full incoming field at each node/label combination.
    std::vector<float> H((size_t)N * K);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < K; ++k) {
            H[(size_t)i * K + k] = (float)node_potential(i, k);
        }
    }

    int iters = 0;
    float max_delta = 0.0f;
    double sum_delta = 0.0;
    double mean_delta = 0.0;

    for (int iter = 0; iter < max_iter; ++iter) {
        max_delta = 0.0f;
        sum_delta = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            std::vector<float> tmp_h(K);
            std::vector<float> tmp_msg(K);
            float thread_max_delta = 0.0f;
            double thread_sum_delta = 0.0;

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int u = 0; u < N; ++u) {
                for (int edge = adj_ptr[u]; edge < adj_ptr[u + 1]; ++edge) {
                    const int rev = rev_idx[edge];

                    const size_t u_offset = (size_t)u * K;
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

                    // Compute next-iteration messages into msg_new.
                    for (int k = 0; k < K; ++k) {
                        const float oldv = msg[edge_offset + k];
                        const float newv = one_minus_damping * oldv + damping * tmp_msg[k];
                        msg_new[edge_offset + k] = newv;

                        const float d = std::fabs(newv - oldv);
                        if (d > thread_max_delta) thread_max_delta = d;
                        thread_sum_delta += d;
                    }
                }
            }

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if (thread_max_delta > max_delta) max_delta = thread_max_delta;
                sum_delta += thread_sum_delta;
            }
        }

        // Rebuild H from scratch from msg_new (serial)
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

        if (verbose) {
            mean_delta = sum_delta / ((double)E * (double)K);
            Rprintf("LBP: Iteration %d, max_delta = %.6g, mean_delta = %.6g\n",
                    iter + 1, max_delta, mean_delta);
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

    mean_delta = sum_delta / ((double)E * (double)K);

    return List::create(
        _["marginals"] = marginals,
        _["iterations"] = iters,
        _["max_delta"] = (double)max_delta,
        _["mean_delta"] = mean_delta
    );
}
