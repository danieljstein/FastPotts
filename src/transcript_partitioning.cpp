// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

static inline std::uint64_t basin_pair_key(const int a0, const int b0) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a0)) << 32) |
        static_cast<std::uint32_t>(b0);
}

static inline std::uint64_t node_pair_key(const int a0, const int b0) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a0)) << 32) |
        static_cast<std::uint32_t>(b0);
}

static inline double safe_prob(const double x, const double eps) {
    return x < eps ? eps : x;
}

static double js_rows(
    const NumericMatrix& posterior,
    const int i0,
    const int j0,
    const double eps
) {
    const int K = posterior.ncol();
    double p_sum = 0.0;
    double q_sum = 0.0;

    for (int k = 0; k < K; ++k) {
        p_sum += safe_prob(posterior(i0, k), eps);
        q_sum += safe_prob(posterior(j0, k), eps);
    }

    double js = 0.0;
    for (int k = 0; k < K; ++k) {
        const double p = safe_prob(posterior(i0, k), eps) / p_sum;
        const double q = safe_prob(posterior(j0, k), eps) / q_sum;
        const double m = 0.5 * (p + q);
        js += 0.5 * p * std::log(p / m) + 0.5 * q * std::log(q / m);
    }
    return js;
}

static int find_root(std::vector<int>& parent, int x) {
    int root = x;
    while (parent[root] != root) {
        root = parent[root];
    }
    while (parent[x] != x) {
        const int next = parent[x];
        parent[x] = root;
        x = next;
    }
    return root;
}

struct EdgeStats {
    int from;
    int to;
    double distance;

    EdgeStats() : from(0), to(0), distance(std::numeric_limits<double>::infinity()) {}
    EdgeStats(const int from_, const int to_, const double distance_) :
        from(from_), to(to_), distance(distance_) {}
};

static inline std::uint64_t type_basin_key(const int type0, const int basin0) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(type0)) << 32) |
        static_cast<std::uint32_t>(basin0);
}

struct BasinMeasureStats {
    int type;
    int basin;
    int n_simplex;
    double measure;

    BasinMeasureStats() : type(0), basin(0), n_simplex(0), measure(0.0) {}
    BasinMeasureStats(const int type_, const int basin_) :
        type(type_), basin(basin_), n_simplex(0), measure(0.0) {}
};

//' Compact directed KNN results into undirected graph edges
//'
//' Internal C++ helper for `partition_transcripts_watershed()`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
DataFrame compact_knn_edges_cpp(
    const IntegerMatrix& nn_idx,
    const NumericMatrix& nn_dist,
    const IntegerVector& query_index,
    const int n_nodes,
    const double max_distance
) {
    const int n_query = nn_idx.nrow();
    const int k = nn_idx.ncol();
    if (nn_dist.nrow() != n_query || nn_dist.ncol() != k) {
        stop("nn_idx and nn_dist must have the same dimensions.");
    }
    if (query_index.size() != n_query) {
        stop("query_index must have one entry per nearest-neighbor query row.");
    }

    std::unordered_map<std::uint64_t, EdgeStats> edges;
    edges.reserve(static_cast<std::size_t>(n_query) * static_cast<std::size_t>(k));

    for (int i = 0; i < n_query; ++i) {
        const int q = query_index[i] - 1;
        if (q < 0 || q >= n_nodes) {
            stop("query_index contains an out-of-range node index.");
        }
        for (int j = 0; j < k; ++j) {
            const int nbr = nn_idx(i, j) - 1;
            const double dist = nn_dist(i, j);
            if (nbr < 0 || nbr >= n_nodes || nbr == q) {
                continue;
            }
            if (!R_finite(dist) || dist > max_distance) {
                continue;
            }

            int a = q;
            int b = nbr;
            if (a > b) {
                std::swap(a, b);
            }
            const std::uint64_t key = node_pair_key(a, b);
            auto it = edges.find(key);
            if (it == edges.end()) {
                edges.emplace(key, EdgeStats(a, b, dist));
            } else if (dist < it->second.distance) {
                it->second.distance = dist;
            }
        }
    }

    std::vector<EdgeStats> rows;
    rows.reserve(edges.size());
    for (const auto& kv : edges) {
        rows.push_back(kv.second);
    }
    std::sort(
        rows.begin(),
        rows.end(),
        [](const EdgeStats& lhs, const EdgeStats& rhs) {
            if (lhs.from != rhs.from) return lhs.from < rhs.from;
            return lhs.to < rhs.to;
        }
    );

    const int E = rows.size();
    IntegerVector from(E);
    IntegerVector to(E);
    NumericVector distance(E);
    for (int e = 0; e < E; ++e) {
        from[e] = rows[e].from + 1;
        to[e] = rows[e].to + 1;
        distance[e] = rows[e].distance;
    }

    return DataFrame::create(
        _["from"] = from,
        _["to"] = to,
        _["distance"] = distance
    );
}

//' Compact an edge list into unique undirected graph edges
//'
//' Internal C++ helper for `partition_transcripts_watershed()`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
DataFrame compact_undirected_edges_cpp(
    const IntegerVector& from,
    const IntegerVector& to,
    const NumericVector& distance,
    const int n_nodes,
    const double max_distance
) {
    const R_xlen_t E_in = from.size();
    if (to.size() != E_in || distance.size() != E_in) {
        stop("from, to, and distance must have the same length.");
    }

    std::unordered_map<std::uint64_t, EdgeStats> edges;
    edges.reserve(static_cast<std::size_t>(E_in));

    for (R_xlen_t e = 0; e < E_in; ++e) {
        int a = from[e] - 1;
        int b = to[e] - 1;
        const double dist = distance[e];
        if (a < 0 || a >= n_nodes || b < 0 || b >= n_nodes || a == b) {
            continue;
        }
        if (!R_finite(dist) || dist > max_distance) {
            continue;
        }
        if (a > b) {
            std::swap(a, b);
        }
        const std::uint64_t key = node_pair_key(a, b);
        auto it = edges.find(key);
        if (it == edges.end()) {
            edges.emplace(key, EdgeStats(a, b, dist));
        } else if (dist < it->second.distance) {
            it->second.distance = dist;
        }
    }

    std::vector<EdgeStats> rows;
    rows.reserve(edges.size());
    for (const auto& kv : edges) {
        rows.push_back(kv.second);
    }
    std::sort(
        rows.begin(),
        rows.end(),
        [](const EdgeStats& lhs, const EdgeStats& rhs) {
            if (lhs.from != rhs.from) return lhs.from < rhs.from;
            return lhs.to < rhs.to;
        }
    );

    const int E = rows.size();
    IntegerVector out_from(E);
    IntegerVector out_to(E);
    NumericVector out_distance(E);
    for (int e = 0; e < E; ++e) {
        out_from[e] = rows[e].from + 1;
        out_to[e] = rows[e].to + 1;
        out_distance[e] = rows[e].distance;
    }

    return DataFrame::create(
        _["from"] = out_from,
        _["to"] = out_to,
        _["distance"] = out_distance
    );
}

//' Aggregate basis-watershed simplex area or volume by basin
//'
//' Internal C++ helper for `basis_watershed_gene_count_data()`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List basis_watershed_domain_measure_cpp(
    const IntegerMatrix& domain_basis_id,
    const NumericVector& simplex_volume,
    const List& active_start,
    const List& type_density,
    const List& basin
) {
    const int n_simplex = domain_basis_id.nrow();
    const int n_vertex = domain_basis_id.ncol();
    const int K = active_start.size();
    if (type_density.size() != K || basin.size() != K) {
        stop("active_start, type_density, and basin must have the same length.");
    }
    if (simplex_volume.size() != 1 && simplex_volume.size() != n_simplex) {
        stop("simplex_volume must have length 1 or one value per simplex.");
    }

    std::vector<LogicalVector> active(K);
    std::vector<NumericVector> density(K);
    std::vector<IntegerVector> basin_vec(K);
    for (int k = 0; k < K; ++k) {
        active[k] = as<LogicalVector>(active_start[k]);
        density[k] = as<NumericVector>(type_density[k]);
        basin_vec[k] = as<IntegerVector>(basin[k]);
        if (density[k].size() != active[k].size() || basin_vec[k].size() != active[k].size()) {
            stop("Each active_start, type_density, and basin vector must have the same length.");
        }
    }

    std::unordered_map<std::uint64_t, BasinMeasureStats> stats;
    stats.reserve(static_cast<std::size_t>(K) * 1024);
    std::vector<int> active_count(K);
    std::vector<int> best_basis0(K);
    std::vector<double> best_density(K);

    double domain_measure_total = 0.0;
    double active_domain_measure = 0.0;
    int n_active_domain_simplexes = 0;

    for (int i = 0; i < n_simplex; ++i) {
        const double vol = simplex_volume.size() == 1 ? simplex_volume[0] : simplex_volume[i];
        if (!R_finite(vol)) {
            stop("simplex_volume must contain finite values.");
        }
        domain_measure_total += vol;

        int total_active_vertices = 0;
        std::fill(active_count.begin(), active_count.end(), 0);
        std::fill(best_basis0.begin(), best_basis0.end(), -1);
        std::fill(best_density.begin(), best_density.end(), R_NegInf);

        for (int v = 0; v < n_vertex; ++v) {
            const int b0 = domain_basis_id(i, v) - 1;
            if (b0 < 0) {
                continue;
            }
            for (int k = 0; k < K; ++k) {
                if (b0 >= active[k].size()) {
                    stop("domain_basis_id contains a basis index outside active_start.");
                }
                const int is_active = active[k][b0];
                if (is_active == TRUE) {
                    active_count[k] += 1;
                    total_active_vertices += 1;
                    const double dens = density[k][b0];
                    if (R_finite(dens) && dens > best_density[k]) {
                        best_density[k] = dens;
                        best_basis0[k] = b0;
                    }
                }
            }
        }

        if (total_active_vertices == 0) {
            continue;
        }
        n_active_domain_simplexes += 1;
        active_domain_measure += vol;

        for (int k = 0; k < K; ++k) {
            if (active_count[k] == 0 || best_basis0[k] < 0) {
                continue;
            }
            const int basin_id = basin_vec[k][best_basis0[k]];
            if (basin_id == NA_INTEGER) {
                continue;
            }
            const double weight = static_cast<double>(active_count[k]) /
                static_cast<double>(total_active_vertices);
            const std::uint64_t key = type_basin_key(k, basin_id - 1);
            auto it = stats.find(key);
            if (it == stats.end()) {
                it = stats.emplace(key, BasinMeasureStats(k + 1, basin_id)).first;
            }
            it->second.n_simplex += 1;
            it->second.measure += vol * weight;
        }
    }

    std::vector<BasinMeasureStats> rows;
    rows.reserve(stats.size());
    for (const auto& kv : stats) {
        rows.push_back(kv.second);
    }
    std::sort(
        rows.begin(),
        rows.end(),
        [](const BasinMeasureStats& lhs, const BasinMeasureStats& rhs) {
            if (lhs.type != rhs.type) return lhs.type < rhs.type;
            return lhs.basin < rhs.basin;
        }
    );

    const int n_row = rows.size();
    IntegerVector type_index(n_row);
    IntegerVector basin_id(n_row);
    IntegerVector n_domain_simplexes(n_row);
    NumericVector basis_measure(n_row);
    for (int r = 0; r < n_row; ++r) {
        type_index[r] = rows[r].type;
        basin_id[r] = rows[r].basin;
        n_domain_simplexes[r] = rows[r].n_simplex;
        basis_measure[r] = rows[r].measure;
    }

    DataFrame measure = DataFrame::create(
        _["type_index"] = type_index,
        _["basin"] = basin_id,
        _["n_domain_simplexes"] = n_domain_simplexes,
        _["basis_measure"] = basis_measure
    );
    DataFrame domain_metadata = DataFrame::create(
        _["n_domain_simplexes_total"] = n_simplex,
        _["domain_measure_total"] = domain_measure_total,
        _["n_active_domain_simplexes"] = n_active_domain_simplexes,
        _["active_domain_measure"] = active_domain_measure
    );

    return List::create(
        _["measure"] = measure,
        _["domain_metadata"] = domain_metadata
    );
}

//' Estimate graph transcript density
//'
//' Internal C++ helper for `partition_transcripts_watershed()`.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector estimate_graph_density_cpp(
    const IntegerVector& from,
    const IntegerVector& to,
    const NumericVector& distance,
    const NumericMatrix& posterior,
    const double bandwidth,
    const int mode
) {
    const int n = posterior.nrow();
    const int K = posterior.ncol();
    const R_xlen_t E = from.size();

    if (to.size() != E || distance.size() != E) {
        stop("from, to, and distance must have the same length.");
    }
    if (!R_finite(bandwidth) || bandwidth <= 0.0) {
        stop("bandwidth must be positive and finite.");
    }

    NumericVector density(n);

    if (mode == 1) {
        std::fill(density.begin(), density.end(), 1.0);
        for (R_xlen_t e = 0; e < E; ++e) {
            const int a = from[e] - 1;
            const int b = to[e] - 1;
            const double r = distance[e] / bandwidth;
            const double kernel = std::exp(-0.5 * r * r);
            density[a] += kernel;
            density[b] += kernel;
        }
        return density;
    }

    NumericMatrix type_density(clone(posterior));
    for (R_xlen_t e = 0; e < E; ++e) {
        const int a = from[e] - 1;
        const int b = to[e] - 1;
        const double r = distance[e] / bandwidth;
        const double kernel = std::exp(-0.5 * r * r);

        for (int k = 0; k < K; ++k) {
            type_density(a, k) += kernel * posterior(b, k);
            type_density(b, k) += kernel * posterior(a, k);
        }
    }

    for (int i = 0; i < n; ++i) {
        double value = 0.0;
        for (int k = 0; k < K; ++k) {
            value += posterior(i, k) * type_density(i, k);
        }
        density[i] = value;
    }

    return density;
}

//' Compute edge-level posterior Jensen-Shannon divergence
//'
//' Internal C++ helper for `partition_transcripts_watershed()`.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector posterior_js_divergence_edges_cpp(
    const NumericMatrix& posterior,
    const IntegerVector& from,
    const IntegerVector& to,
    const double eps = 1e-12
) {
    const R_xlen_t E = from.size();
    if (to.size() != E) {
        stop("from and to must have the same length.");
    }

    NumericVector out(E);
    for (R_xlen_t e = 0; e < E; ++e) {
        out[e] = js_rows(posterior, from[e] - 1, to[e] - 1, eps);
    }
    return out;
}

//' Assign graph nodes to density modes by ascent
//'
//' Internal C++ helper for `partition_transcripts_watershed()`.
//'
//' @keywords internal
// [[Rcpp::export]]
List density_ascent_partition_cpp(
    const IntegerVector& from,
    const IntegerVector& to,
    const NumericVector& distance,
    const NumericVector& density,
    const NumericVector& posterior_js,
    const double distance_weight,
    const double posterior_weight
) {
    const int n = density.size();
    const R_xlen_t E = from.size();
    if (to.size() != E || distance.size() != E || posterior_js.size() != E) {
        stop("from, to, distance, and posterior_js must have the same length.");
    }

    std::vector<int> parent(n);
    std::vector<double> best_score(n, -std::numeric_limits<double>::infinity());
    for (int i = 0; i < n; ++i) {
        parent[i] = i;
    }

    auto update_parent = [&](const int i, const int j, const double dist, const double js) {
        if (density[j] <= density[i]) {
            return;
        }
        const double score = (density[j] - density[i]) -
            distance_weight * dist -
            posterior_weight * js;
        if (score > best_score[i]) {
            best_score[i] = score;
            parent[i] = j;
        }
    };

    for (R_xlen_t e = 0; e < E; ++e) {
        const int a = from[e] - 1;
        const int b = to[e] - 1;
        update_parent(a, b, distance[e], posterior_js[e]);
        update_parent(b, a, distance[e], posterior_js[e]);
    }

    std::vector<int> roots(n);
    std::vector<int> unique_roots;
    unique_roots.reserve(n);
    for (int i = 0; i < n; ++i) {
        roots[i] = find_root(parent, i);
        unique_roots.push_back(roots[i]);
    }

    std::sort(unique_roots.begin(), unique_roots.end());
    unique_roots.erase(std::unique(unique_roots.begin(), unique_roots.end()), unique_roots.end());

    IntegerVector parent_out(n);
    IntegerVector basin(n);
    for (int i = 0; i < n; ++i) {
        parent_out[i] = parent[i] + 1;
        basin[i] = static_cast<int>(
            std::lower_bound(unique_roots.begin(), unique_roots.end(), roots[i]) - unique_roots.begin()
        ) + 1;
    }

    return List::create(
        _["parent"] = parent_out,
        _["basin"] = basin
    );
}

//' Assign active graph nodes to density modes by ascent
//'
//' Internal C++ helper for basis-level watershed partitioning. Parents are
//' computed over all graph nodes, but basin labels are assigned only for
//' active starts. Inactive nodes can therefore be used as pass-through nodes
//' along ascent paths without producing their own output basins.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List density_ascent_active_partition_cpp(
    const IntegerVector& from,
    const IntegerVector& to,
    const NumericVector& distance,
    const NumericVector& density,
    const LogicalVector& active_start,
    const double distance_weight
) {
    const int n = density.size();
    const R_xlen_t E = from.size();
    if (to.size() != E || distance.size() != E) {
        stop("from, to, and distance must have the same length.");
    }
    if (active_start.size() != n) {
        stop("active_start must have one value per graph node.");
    }

    std::vector<int> parent(n);
    std::vector<double> best_score(n, -std::numeric_limits<double>::infinity());
    for (int i = 0; i < n; ++i) {
        parent[i] = i;
    }

    auto update_parent = [&](const int i, const int j, const double dist) {
        if (density[j] <= density[i]) {
            return;
        }
        const double score = (density[j] - density[i]) - distance_weight * dist;
        if (score > best_score[i]) {
            best_score[i] = score;
            parent[i] = j;
        }
    };

    for (R_xlen_t e = 0; e < E; ++e) {
        const int a = from[e] - 1;
        const int b = to[e] - 1;
        update_parent(a, b, distance[e]);
        update_parent(b, a, distance[e]);
    }

    std::vector<int> root(n);
    for (int i = 0; i < n; ++i) {
        root[i] = find_root(parent, i);
    }

    std::vector<int> unique_active_roots;
    unique_active_roots.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (active_start[i] == TRUE) {
            unique_active_roots.push_back(root[i]);
        }
    }
    std::sort(unique_active_roots.begin(), unique_active_roots.end());
    unique_active_roots.erase(
        std::unique(unique_active_roots.begin(), unique_active_roots.end()),
        unique_active_roots.end()
    );

    IntegerVector parent_out(n);
    IntegerVector root_out(n);
    IntegerVector basin(n, NA_INTEGER);
    for (int i = 0; i < n; ++i) {
        parent_out[i] = parent[i] + 1;
        root_out[i] = root[i] + 1;
        if (active_start[i] == TRUE) {
            basin[i] = static_cast<int>(
                std::lower_bound(unique_active_roots.begin(), unique_active_roots.end(), root[i]) -
                unique_active_roots.begin()
            ) + 1;
        }
    }

    IntegerVector mode_basis(unique_active_roots.size());
    for (std::size_t i = 0; i < unique_active_roots.size(); ++i) {
        mode_basis[i] = unique_active_roots[i] + 1;
    }

    return List::create(
        _["parent"] = parent_out,
        _["root"] = root_out,
        _["basin"] = basin,
        _["mode_basis"] = mode_basis
    );
}

struct BasinPairStats {
    int a;
    int b;
    int edge_count;
    double saddle_density;
    double distance_sum;
    double min_distance;
    double edge_js_sum;

    BasinPairStats() :
        a(0),
        b(0),
        edge_count(0),
        saddle_density(-std::numeric_limits<double>::infinity()),
        distance_sum(0.0),
        min_distance(std::numeric_limits<double>::infinity()),
        edge_js_sum(0.0) {}
};

//' Build basin adjacency and saddle diagnostics
//'
//' Internal C++ helper for `partition_transcripts_watershed()`.
//'
//' @keywords internal
// [[Rcpp::export]]
DataFrame build_basin_adjacency_cpp(
    const IntegerVector& from,
    const IntegerVector& to,
    const NumericVector& distance,
    const IntegerVector& basin,
    const NumericVector& density,
    const NumericMatrix& posterior,
    const NumericVector& posterior_js,
    const double eps = 1e-12
) {
    const int n = basin.size();
    const int K = posterior.ncol();
    const R_xlen_t E = from.size();
    if (to.size() != E || distance.size() != E || posterior_js.size() != E) {
        stop("from, to, distance, and posterior_js must have the same length.");
    }
    if (density.size() != n || posterior.nrow() != n) {
        stop("basin, density, and posterior dimensions are inconsistent.");
    }

    int B = 0;
    for (int i = 0; i < n; ++i) {
        if (basin[i] > B) {
            B = basin[i];
        }
    }

    std::vector<double> mode_density(B, -std::numeric_limits<double>::infinity());
    std::vector<double> posterior_sum(static_cast<std::size_t>(B) * K, 0.0);
    std::vector<double> posterior_total(B, 0.0);

    for (int i = 0; i < n; ++i) {
        const int b = basin[i] - 1;
        if (density[i] > mode_density[b]) {
            mode_density[b] = density[i];
        }
        for (int k = 0; k < K; ++k) {
            const double value = posterior(i, k);
            posterior_sum[static_cast<std::size_t>(b) * K + k] += value;
            posterior_total[b] += value;
        }
    }

    std::unordered_map<std::uint64_t, BasinPairStats> stats;
    stats.reserve(static_cast<std::size_t>(E / 4 + 1));

    for (R_xlen_t e = 0; e < E; ++e) {
        int a = basin[from[e] - 1] - 1;
        int b = basin[to[e] - 1] - 1;
        if (a == b) {
            continue;
        }
        if (a > b) {
            std::swap(a, b);
        }

        const std::uint64_t key = basin_pair_key(a, b);
        auto it = stats.find(key);
        if (it == stats.end()) {
            BasinPairStats init;
            init.a = a;
            init.b = b;
            it = stats.emplace(key, init).first;
        }

        BasinPairStats& s = it->second;
        const double saddle = std::min(density[from[e] - 1], density[to[e] - 1]);
        s.edge_count += 1;
        if (saddle > s.saddle_density) {
            s.saddle_density = saddle;
        }
        s.distance_sum += distance[e];
        if (distance[e] < s.min_distance) {
            s.min_distance = distance[e];
        }
        s.edge_js_sum += posterior_js[e];
    }

    std::vector<BasinPairStats> rows;
    rows.reserve(stats.size());
    for (const auto& kv : stats) {
        rows.push_back(kv.second);
    }
    std::sort(
        rows.begin(),
        rows.end(),
        [](const BasinPairStats& lhs, const BasinPairStats& rhs) {
            if (lhs.a != rhs.a) return lhs.a < rhs.a;
            return lhs.b < rhs.b;
        }
    );

    const int R = rows.size();
    IntegerVector basin_a(R);
    IntegerVector basin_b(R);
    IntegerVector edge_count(R);
    NumericVector saddle_density(R);
    NumericVector saddle_ratio(R);
    NumericVector mean_distance(R);
    NumericVector min_distance(R);
    NumericVector mean_edge_js(R);
    NumericVector basin_js(R);

    for (int r = 0; r < R; ++r) {
        const BasinPairStats& s = rows[r];
        basin_a[r] = s.a + 1;
        basin_b[r] = s.b + 1;
        edge_count[r] = s.edge_count;
        saddle_density[r] = s.saddle_density;
        saddle_ratio[r] = s.saddle_density / std::min(mode_density[s.a], mode_density[s.b]);
        mean_distance[r] = s.distance_sum / s.edge_count;
        min_distance[r] = s.min_distance;
        mean_edge_js[r] = s.edge_js_sum / s.edge_count;

        double js = 0.0;
        for (int k = 0; k < K; ++k) {
            const double p = safe_prob(
                posterior_sum[static_cast<std::size_t>(s.a) * K + k] / posterior_total[s.a],
                eps
            );
            const double q = safe_prob(
                posterior_sum[static_cast<std::size_t>(s.b) * K + k] / posterior_total[s.b],
                eps
            );
            const double m = 0.5 * (p + q);
            js += 0.5 * p * std::log(p / m) + 0.5 * q * std::log(q / m);
        }
        basin_js[r] = js;
    }

    return DataFrame::create(
        _["basin_a"] = basin_a,
        _["basin_b"] = basin_b,
        _["edge_count"] = edge_count,
        _["saddle_density"] = saddle_density,
        _["saddle_ratio"] = saddle_ratio,
        _["mean_distance"] = mean_distance,
        _["min_distance"] = min_distance,
        _["mean_edge_js"] = mean_edge_js,
        _["basin_js"] = basin_js
    );
}
