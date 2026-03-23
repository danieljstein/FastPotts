// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

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
    const int E = adj_idx.size();          // number of edges

    if (adj_ptr.size() != N + 1) stop("adj_ptr must have length N+1");
    if (rev_idx.size() != E || edge_weights.size() != E) stop("edge vectors must match");
    if (damping <= 0.0 || damping > 1.0) stop("damping must be in (0,1]");

    const float one_minus_damping = 1.0f - damping;

    std::vector<float> w(E);    // edge weights in float: log of same-label vs. different-label assignment
    for (int e = 0; e < E; ++e) {
        w[e] = (float) edge_weights[e];
    }

    // log-messages on directed edges
    std::vector<float> msg((size_t)E * K, 0.0f);
    std::vector<float> msg_new;
    if (synchronous) {
        msg_new.assign((size_t)E * K, 0.0f);
    }

    // H_i(k) = log phi_i(k) + sum_{l in N(i)} log m_{l->i}(k)
    std::vector<float> H((size_t)N * K);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < K; ++k) {
            H[(size_t)i * K + k] = (float)node_potential(i, k);
        }
    }

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

                // cavity field h_{u->v}(k) = H_u(k) - log m_{v->u}(k)
                for (int k = 0; k < K; ++k) {
                    tmp_h[k] = H[u_offset + k] - msg[rev_offset + k];
                }

                // log S where S = sum_a exp(h[a])
                const float logS = logsumexp_float(tmp_h.data(), K);
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

                // normalize message in log-domain
                const float lse = logsumexp_float(tmp_msg.data(), K);
                for (int k = 0; k < K; ++k) {
                    tmp_msg[k] -= lse;
                }

                if (synchronous) {
                    for (int k = 0; k < K; ++k) {
                        const float oldv = msg[edge_offset + k];
                        const float newv = one_minus_damping * oldv + damping * tmp_msg[k];
                        msg_new[edge_offset + k] = newv;

                        const float d = std::fabs(newv - oldv);
                        if (d > max_delta) max_delta = d;
                    }
                } else {
                    // in-place asynchronous update
                    for (int k = 0; k < K; ++k) {
                        const float oldv = msg[edge_offset + k];
                        const float newv = one_minus_damping * oldv + damping * tmp_msg[k];
                        msg[edge_offset + k] = newv;

                        // message u->v contributes to H[v]
                        H[v_offset + k] += (newv - oldv);

                        const float d = std::fabs(newv - oldv);
                        if (d > max_delta) max_delta = d;
                    }
                }
            }
        }

        if (synchronous) {
            // rebuild H from scratch from msg_new
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
