#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

using namespace std;

namespace hollan {

struct Problem {
    int rows = 0, cols = 0, n = 0;
    vector<int> rank_at;
    vector<int> pos;
    vector<vector<int>> nbr;

    explicit Problem(const vector<vector<int>>& matrix) {
        rows = static_cast<int>(matrix.size());
        cols = rows ? static_cast<int>(matrix[0].size()) : 0;
        n = rows * cols;
        rank_at.assign(n, 0);
        pos.assign(n, 0);
        vector<int> cells(n);
        iota(cells.begin(), cells.end(), 0);
        sort(cells.begin(), cells.end(), [&](int a, int b) {
            int va = matrix[a / cols][a % cols];
            int vb = matrix[b / cols][b % cols];
            return va != vb ? va < vb : a < b;
        });
        for (int k = 0; k < n; ++k) {
            pos[k] = cells[k];
            rank_at[cells[k]] = k;
        }
        nbr.assign(n, {});
        static const int dr[4] = {-1, 1, 0, 0};
        static const int dc[4] = {0, 0, -1, 1};
        for (int v = 0; v < n; ++v) {
            int r = v / cols, c = v % cols;
            for (int d = 0; d < 4; ++d) {
                int nr = r + dr[d], nc = c + dc[d];
                if (0 <= nr && nr < rows && 0 <= nc && nc < cols)
                    nbr[v].push_back(nr * cols + nc);
            }
        }
    }
};

struct Chains {
    vector<int> parent, size, start, finish;

    explicit Chains(int n = 0) : parent(n), size(n, 1), start(n), finish(n) {
        iota(parent.begin(), parent.end(), 0);
        iota(start.begin(), start.end(), 0);
        iota(finish.begin(), finish.end(), 0);
    }

    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    bool add_edge(int from, int to) {
        int a = find(from), b = find(to);
        if (a == b) return true;
        int new_start = start[a], new_finish = finish[b];
        if (size[a] < size[b]) swap(a, b);
        parent[b] = a;
        size[a] += size[b];
        start[a] = new_start;
        finish[a] = new_finish;
        return false;
    }
};

struct Weights {
    double close = 10000.0;
    double reserved = 300.0;
    double ready = 20.0;
    double expose = 2.0;
    double new_frontier = 1.0;
    double selected_neighbors = 0.3;
    double future_rank = 0.0;
    double finish_delay = 0.0;
    double noise = 0.0;
    int lookahead = 12;
};

struct Candidate {
    vector<int> order;
    int cycles = -1;
    int worst = numeric_limits<int>::max();
};

// Populated only by offline benchmark tools. Passing nullptr has no effect on
// the submitted solver's decisions or output.
struct SearchStats {
    int fallback_swaps = -1;
    int deterministic_swaps = -1;
    int post_local_swaps = -1;
    int post_beam_swaps = -1;
    int post_exact_swaps = -1;
    int final_swaps = -1;
    int deterministic_wins = 0;
    int local_wins = 0;
    int beam_wins = 0;
    int exact_wins = 0;
    int random_beam_wins = 0;
    int random_repair_wins = 0;
    int random_restart_wins = 0;
    int destroy_repair_wins = 0;
    int iterations = 0;
    double deterministic_ms = 0.0;
    double local_ms = 0.0;
    double beam_ms = 0.0;
    double exact_ms = 0.0;
    double last_improvement_ms = 0.0;
};

// Guaranteed connected-prefix ordering. It fixes every rank whose original
// cell is already on the frontier and otherwise displaces a late rank.
static vector<int> baseline_order(const Problem& p) {
    vector<int> order;
    order.reserve(p.n);
    if (p.n == 0) return order;
    vector<unsigned char> used(p.n, 0), in_frontier(p.n, 0);
    vector<int> frontier;
    auto add_cell = [&](int v) {
        used[v] = 1;
        order.push_back(v);
        for (int u : p.nbr[v]) {
            if (!used[u] && !in_frontier[u]) {
                in_frontier[u] = 1;
                frontier.push_back(u);
            }
        }
    };
    add_cell(p.pos[0]);
    for (int k = 1; k < p.n; ++k) {
        int chosen = -1;
        if (!used[p.pos[k]] && in_frontier[p.pos[k]]) {
            chosen = p.pos[k];
        } else {
            for (int v : frontier) {
                if (!used[v] && (chosen < 0 || p.rank_at[v] > p.rank_at[chosen]))
                    chosen = v;
            }
        }
        in_frontier[chosen] = 0;
        add_cell(chosen);
    }
    return order;
}

static int cycle_count(const Problem& p, const vector<int>& order) {
    vector<unsigned char> seen(p.n, 0);
    int cycles = 0;
    for (int s = 0; s < p.n; ++s) if (!seen[s]) {
        ++cycles;
        for (int v = s; !seen[v]; v = p.rank_at[order[v]]) seen[v] = 1;
    }
    return cycles;
}

static int worst_distance(const Problem& p, const vector<int>& order) {
    vector<int> target_rank(p.n), dist(p.n, 0);
    for (int k = 0; k < p.n; ++k) target_rank[order[k]] = k;
    int answer = 0;
    for (int k = 1; k < p.n; ++k) {
        int v = order[k], best = -1;
        for (int u : p.nbr[v]) if (target_rank[u] < k)
            best = max(best, dist[u]);
        dist[v] = best + 1;
        answer = max(answer, dist[v]);
    }
    return answer;
}

static Candidate evaluate(const Problem& p, vector<int> order) {
    Candidate c;
    c.cycles = cycle_count(p, order);
    c.worst = worst_distance(p, order);
    c.order = std::move(order);
    return c;
}

static bool better(const Candidate& a, const Candidate& b) {
    return a.cycles > b.cycles || (a.cycles == b.cycles && a.worst < b.worst);
}

class Builder {
public:
    const Problem& p;
    mt19937_64& rng;
    Weights w;
    vector<unsigned char> used, in_frontier;
    vector<int> frontier, order;
    Chains chains;
    int cycles = 0;

    Builder(const Problem& problem, mt19937_64& generator, const Weights& weights)
        : p(problem), rng(generator), w(weights), used(p.n, 0),
          in_frontier(p.n, 0), chains(p.n) {
        order.reserve(p.n);
        frontier.reserve(p.n);
    }

    void add_cell(int v) {
        int k = static_cast<int>(order.size());
        used[v] = 1;
        in_frontier[v] = 0;
        order.push_back(v);
        cycles += chains.add_edge(k, p.rank_at[v]);
        for (int u : p.nbr[v]) {
            if (!used[u] && !in_frontier[u]) {
                in_frontier[u] = 1;
                frontier.push_back(u);
            }
        }
    }

    bool replay_prefix(const vector<int>& prefix, int length) {
        for (int k = 0; k < length; ++k) {
            int v = prefix[k];
            if (used[v]) return false;
            if (k && !in_frontier[v]) return false;
            add_cell(v);
        }
        return true;
    }

    double score(int v, int k) {
        int label = p.rank_at[v];
        int rk = chains.find(k), ra = chains.find(label);
        bool closes = rk == ra;
        int selected_nbr = 0, fresh = 0;
        double expose_score = 0.0;
        for (int u : p.nbr[v]) {
            if (used[u]) {
                ++selected_nbr;
            } else if (!in_frontier[u]) {
                ++fresh;
                int d = p.rank_at[u] - k;
                if (1 <= d && d <= w.lookahead)
                    expose_score += 1.0 / d;
            }
        }

        bool is_reserved = false;
        if (!closes && chains.size[ra] > 1) {
            int endpoint = chains.finish[ra];
            is_reserved = endpoint > k;
        }

        int chain_start = chains.start[rk];
        int chain_finish = chains.finish[ra];
        int start_cell = p.pos[chain_start];
        bool ready = !closes && (in_frontier[start_cell] ||
                     find(p.nbr[v].begin(), p.nbr[v].end(), start_cell) != p.nbr[v].end());

        double result = (closes ? w.close : 0.0)
                      - (is_reserved ? w.reserved : 0.0)
                      + (ready ? w.ready : 0.0)
                      + w.expose * expose_score
                      + w.new_frontier * fresh
                      + w.selected_neighbors * selected_nbr
                      + w.future_rank * (double(label - k) / max(1, p.n))
                      + w.finish_delay * (double(chain_finish - k) / max(1, p.n));
        if (w.noise > 0.0) {
            double unit = double(rng() >> 11) * (1.0 / 9007199254740992.0);
            result += w.noise * (2.0 * unit - 1.0);
        }
        return result;
    }

    vector<int> finish() {
        while (static_cast<int>(order.size()) < p.n) {
            int k = static_cast<int>(order.size());
            int chosen = -1;
            double best_score = -numeric_limits<double>::infinity();
            for (int v : frontier) if (!used[v]) {
                double s = score(v, k);
                if (s > best_score || (s == best_score && v < chosen)) {
                    best_score = s;
                    chosen = v;
                }
            }
            add_cell(chosen);
        }
        return std::move(order);
    }
};

static vector<int> greedy_order(const Problem& p, mt19937_64& rng,
                                const Weights& w, int root,
                                const vector<int>* prefix = nullptr,
                                int prefix_length = 0) {
    Builder b(p, rng, w);
    if (prefix && prefix_length > 0) {
        b.replay_prefix(*prefix, prefix_length);
    } else {
        b.add_cell(root);
    }
    return b.finish();
}

struct RepairWeights {
    double expose_displaced = 100.0;
    double expose_upcoming = 5.0;
    double new_frontier = 1.0;
    double selected_neighbors = 0.0;
    double delay = 0.0;
    double noise = 0.0;
    int lookahead = 16;
};

// Start with the identity target (zero swaps). Whenever its next cell is not
// connected, transpose that entry with a frontier entry. This keeps an exact,
// explicit repair certificate and tends to close the resulting 2-cycles when
// the displaced cell is exposed before its new position is reached.
static vector<int> repair_order(const Problem& p, mt19937_64& rng,
                                const RepairWeights& w) {
    vector<int> seq(p.n), where(p.n), order;
    iota(seq.begin(), seq.end(), 0);
    iota(where.begin(), where.end(), 0);
    order.reserve(p.n);
    vector<unsigned char> used(p.n, 0), in_frontier(p.n, 0);
    vector<int> frontier;
    frontier.reserve(p.n);
    vector<double> upcoming_at(p.n, 0.0);
    vector<int> upcoming_stamp(p.n, 0);
    int stamp = 0;

    auto select = [&](int v) {
        used[v] = 1;
        in_frontier[v] = 0;
        order.push_back(v);
        for (int u : p.nbr[v]) if (!used[u] && !in_frontier[u]) {
            in_frontier[u] = 1;
            frontier.push_back(u);
        }
    };

    // Fixing rank zero is always legal and preserves a cycle.
    select(p.pos[0]);
    for (int k = 1; k < p.n; ++k) {
        int current_cell = p.pos[seq[k]];
        if (!in_frontier[current_cell]) {
            int chosen_label = -1;
            double best_score = -numeric_limits<double>::infinity();
            int displaced = seq[k];
            int displaced_cell = p.pos[displaced];
            ++stamp;
            for (int d = 1; d <= w.lookahead && k + d < p.n; ++d) {
                int u = p.pos[seq[k + d]];
                if (!used[u] && !in_frontier[u]) {
                    upcoming_at[u] = 1.0 / d;
                    upcoming_stamp[u] = stamp;
                }
            }
            for (int v : frontier) if (!used[v]) {
                int label = p.rank_at[v];
                int j = where[label];
                if (j <= k) continue;
                int fresh = 0, selected_nbr = 0;
                double upcoming = 0.0;
                bool exposes_displaced = false;
                for (int u : p.nbr[v]) {
                    if (used[u]) ++selected_nbr;
                    else if (!in_frontier[u]) ++fresh;
                    if (u == displaced_cell) exposes_displaced = true;
                    if (upcoming_stamp[u] == stamp) upcoming += upcoming_at[u];
                }
                double score = (exposes_displaced ? w.expose_displaced : 0.0)
                             + w.expose_upcoming * upcoming
                             + w.new_frontier * fresh
                             + w.selected_neighbors * selected_nbr
                             + w.delay * double(j - k) / max(1, p.n);
                if (w.noise > 0.0) {
                    double unit = double(rng() >> 11) * (1.0 / 9007199254740992.0);
                    score += w.noise * (2.0 * unit - 1.0);
                }
                if (score > best_score || (score == best_score && label < chosen_label)) {
                    best_score = score;
                    chosen_label = label;
                }
            }
            int j = where[chosen_label];
            int other = seq[j];
            swap(seq[k], seq[j]);
            where[displaced] = j;
            where[other] = k;
            current_cell = p.pos[seq[k]];
        }
        select(current_cell);
    }
    return order;
}

struct RepairBeamState {
    vector<int> seq, where, order, frontier;
    vector<unsigned char> used, in_frontier;
    int cycles = 0;
    double heuristic = 0.0, merit = 0.0;

    explicit RepairBeamState(int n = 0)
        : seq(n), where(n), used(n, 0), in_frontier(n, 0), cycles(n) {
        iota(seq.begin(), seq.end(), 0);
        iota(where.begin(), where.end(), 0);
        order.reserve(n);
        frontier.reserve(n);
    }
};

static void repair_beam_select(const Problem& p, RepairBeamState& s, int v) {
    s.used[v] = 1;
    s.in_frontier[v] = 0;
    s.order.push_back(v);
    for (int u : p.nbr[v]) if (!s.used[u] && !s.in_frontier[u]) {
        s.in_frontier[u] = 1;
        s.frontier.push_back(u);
    }
}

// Beam search over the identity-repair process. A state branches only when the
// next identity assignment is disconnected, so it explores a much narrower
// and more relevant space than the general connected-order beam.
static vector<int> repair_beam_order(const Problem& p, const RepairWeights& w,
                                     int width, int branch,
                                     const chrono::steady_clock::time_point& deadline,
                                     bool accumulate_merit = true) {
    if (p.n == 0) return {};
    RepairBeamState initial(p.n);
    repair_beam_select(p, initial, p.pos[0]);
    vector<RepairBeamState> beam{std::move(initial)};
    vector<double> upcoming_at(p.n, 0.0);
    vector<int> upcoming_stamp(p.n, 0);
    int stamp = 0;
    for (int k = 1; k < p.n; ++k) {
        if (chrono::steady_clock::now() >= deadline) return {};
        vector<RepairBeamState> next;
        next.reserve(width * branch);
        for (const RepairBeamState& state : beam) {
            int current_cell = p.pos[state.seq[k]];
            if (state.in_frontier[current_cell]) {
                RepairBeamState child = state;
                repair_beam_select(p, child, current_cell);
                child.merit = 1000000.0 * child.cycles
                            + (accumulate_merit ? child.heuristic : 0.0);
                next.push_back(std::move(child));
                continue;
            }
            vector<pair<double, int>> choices;
            int displaced = state.seq[k];
            int displaced_cell = p.pos[displaced];
            ++stamp;
            for (int d = 1; d <= w.lookahead && k + d < p.n; ++d) {
                int u = p.pos[state.seq[k + d]];
                if (!state.used[u] && !state.in_frontier[u]) {
                    upcoming_at[u] = 1.0 / d;
                    upcoming_stamp[u] = stamp;
                }
            }
            for (int v : state.frontier) if (!state.used[v]) {
                int label = p.rank_at[v], j = state.where[label];
                if (j <= k) continue;
                int fresh = 0, selected_nbr = 0;
                double upcoming = 0.0;
                bool exposes_displaced = false;
                for (int u : p.nbr[v]) {
                    if (state.used[u]) ++selected_nbr;
                    else if (!state.in_frontier[u]) ++fresh;
                    if (u == displaced_cell) exposes_displaced = true;
                    if (upcoming_stamp[u] == stamp) upcoming += upcoming_at[u];
                }
                double score = (exposes_displaced ? w.expose_displaced : 0.0)
                             + w.expose_upcoming * upcoming
                             + w.new_frontier * fresh
                             + w.selected_neighbors * selected_nbr
                             + w.delay * double(j - k) / max(1, p.n);
                choices.push_back({score, label});
            }
            int keep = min(branch, static_cast<int>(choices.size()));
            partial_sort(choices.begin(), choices.begin() + keep, choices.end(),
                         [](const auto& a, const auto& b) { return a.first > b.first; });
            for (int q = 0; q < keep; ++q) {
                RepairBeamState child = state;
                int label = choices[q].second;
                int j = child.where[label];
                int displaced_label = child.seq[k];
                // Swapping two images splits a cycle iff their source indices
                // are currently in the same permutation cycle.
                int x = k;
                do {
                    x = child.seq[x];
                } while (x != k && x != j);
                bool same_cycle = x == j;
                swap(child.seq[k], child.seq[j]);
                child.where[displaced_label] = j;
                child.where[label] = k;
                child.cycles += same_cycle ? 1 : -1;
                repair_beam_select(p, child, p.pos[child.seq[k]]);
                child.heuristic = accumulate_merit
                    ? child.heuristic + choices[q].first : choices[q].first;
                child.merit = 1000000.0 * child.cycles + child.heuristic;
                next.push_back(std::move(child));
            }
        }
        int keep = min(width, static_cast<int>(next.size()));
        if (keep == 0) return {};
        partial_sort(next.begin(), next.begin() + keep, next.end(),
                     [](const RepairBeamState& a, const RepairBeamState& b) {
                         return a.merit > b.merit;
                     });
        next.resize(keep);
        beam = std::move(next);
    }
    int best = 0;
    for (int i = 1; i < static_cast<int>(beam.size()); ++i)
        if (beam[i].cycles > beam[best].cycles) best = i;
    return std::move(beam[best].order);
}

// Swapping adjacent targets changes the cycle count by exactly one. If the two
// source ranks are in the same permutation cycle it splits that cycle, saving
// one swap. The connectivity test below is the only validity condition needed.
static void adjacent_improve(const Problem& p, Candidate& candidate) {
    bool changed = true;
    vector<int> cycle_id(p.n);
    while (changed) {
        changed = false;
        vector<unsigned char> seen(p.n, 0);
        int id = 0;
        for (int s = 0; s < p.n; ++s) if (!seen[s]) {
            for (int v = s; !seen[v]; v = p.rank_at[candidate.order[v]]) {
                seen[v] = 1;
                cycle_id[v] = id;
            }
            ++id;
        }
        vector<int> position(p.n);
        for (int i = 0; i < p.n; ++i) position[candidate.order[i]] = i;
        for (int i = 1; i + 1 < p.n; ++i) {
            if (cycle_id[i] != cycle_id[i + 1]) continue;
            int later = candidate.order[i + 1];
            bool touches_earlier = false;
            for (int u : p.nbr[later]) if (position[u] < i) {
                touches_earlier = true;
                break;
            }
            if (!touches_earlier) continue;
            swap(candidate.order[i], candidate.order[i + 1]);
            ++candidate.cycles;
            changed = true;
            break;
        }
    }
    candidate.worst = worst_distance(p, candidate.order);
}

static bool valid_target_transposition(const Problem& p, const vector<int>& order,
                                       const vector<int>& position, int i, int j) {
    int x = order[i], y = order[j];
    for (int t = i; t <= j; ++t) {
        int v = (t == i ? y : (t == j ? x : order[t]));
        bool touches = (t == 0);
        for (int u : p.nbr[v]) {
            int q = position[u];
            if (q == i) q = j;
            else if (q == j) q = i;
            if (q < t) {
                touches = true;
                break;
            }
        }
        if (!touches) return false;
    }
    return true;
}

// A transposition of two images in the same permutation cycle splits it and
// saves exactly one swap. Search that mathematically improving neighborhood,
// accepting only moves whose affected connected-prefix interval stays valid.
static void transposition_improve(const Problem& p, Candidate& candidate,
                                  mt19937_64& rng, int failed_limit,
                                  const chrono::steady_clock::time_point& deadline) {
    if (p.n < 3) return;
    int failed = 0;
    vector<int> cycle_id(p.n), position(p.n);
    while (failed < failed_limit && chrono::steady_clock::now() < deadline) {
        vector<unsigned char> seen(p.n, 0);
        int id = 0;
        for (int s = 0; s < p.n; ++s) if (!seen[s]) {
            for (int v = s; !seen[v]; v = p.rank_at[candidate.order[v]]) {
                seen[v] = 1;
                cycle_id[v] = id;
            }
            ++id;
        }
        for (int t = 0; t < p.n; ++t) position[candidate.order[t]] = t;

        bool accepted = false;
        int batch = min(failed_limit - failed, max(64, 8 * p.n));
        for (int attempt = 0; attempt < batch; ++attempt) {
            int i = int(rng() % (p.n - 1));
            int j = i + 1 + int(rng() % (p.n - i - 1));
            ++failed;
            if (cycle_id[i] != cycle_id[j]) continue;
            if (!valid_target_transposition(p, candidate.order, position, i, j)) continue;
            swap(candidate.order[i], candidate.order[j]);
            ++candidate.cycles;
            failed = 0;
            accepted = true;
            break;
        }
        if (!accepted && batch == 0) break;
    }
    candidate.worst = worst_distance(p, candidate.order);
}

static uint64_t matrix_seed(const Problem& p) {
    uint64_t h = 0x9e3779b97f4a7c15ULL;
    for (int x : p.rank_at) {
        h ^= uint64_t(x + 0x9e37) + (h << 6) + (h >> 2);
        h *= 0xbf58476d1ce4e5b9ULL;
    }
    return h;
}

struct BeamNode {
    vector<int> order;
    vector<unsigned char> used, in_frontier;
    vector<int> frontier;
    Chains chains;
    int cycles = 0;
    double merit = 0.0;

    explicit BeamNode(int n = 0)
        : used(n, 0), in_frontier(n, 0), chains(n) {
        order.reserve(n);
        frontier.reserve(n);
    }
};

static void beam_add(const Problem& p, BeamNode& s, int v) {
    int k = static_cast<int>(s.order.size());
    s.used[v] = 1;
    s.in_frontier[v] = 0;
    s.order.push_back(v);
    s.cycles += s.chains.add_edge(k, p.rank_at[v]);
    for (int u : p.nbr[v]) if (!s.used[u] && !s.in_frontier[u]) {
        s.in_frontier[u] = 1;
        s.frontier.push_back(u);
    }
}

static double beam_merit(const Problem& p, BeamNode& s) {
    int k = static_cast<int>(s.order.size());
    int promised = 0, fixed = 0;
    for (int root = 0; root < p.n; ++root) if (s.chains.find(root) == root) {
        int a = s.chains.start[root], z = s.chains.finish[root];
        if (z >= k && s.in_frontier[p.pos[a]]) {
            if (s.chains.size[root] > 1) ++promised;
            else if (a == z) ++fixed;
        }
    }
    int active_frontier = 0;
    for (int v : s.frontier) active_frontier += !s.used[v];
    return 1000000.0 * s.cycles + 2000.0 * promised
         + 30.0 * fixed + active_frontier;
}

static vector<int> beam_order(const Problem& p, int width, int branch,
                              const chrono::steady_clock::time_point& deadline) {
    BeamNode initial(p.n);
    beam_add(p, initial, p.pos[0]);
    vector<BeamNode> beam;
    beam.push_back(std::move(initial));
    for (int k = 1; k < p.n; ++k) {
        if (chrono::steady_clock::now() >= deadline) return {};
        vector<BeamNode> next;
        next.reserve(width * branch);
        for (BeamNode& state : beam) {
            vector<pair<double, int>> choices;
            for (int v : state.frontier) if (!state.used[v]) {
                int label = p.rank_at[v];
                int a = state.chains.find(k), b = state.chains.find(label);
                bool closes = a == b;
                bool reserved = !closes && state.chains.size[b] > 1
                             && state.chains.finish[b] > k;
                int fresh = 0, selected = 0;
                for (int u : p.nbr[v]) {
                    fresh += !state.used[u] && !state.in_frontier[u];
                    selected += state.used[u];
                }
                double score = (closes ? 100000.0 : 0.0)
                             - (reserved ? 5000.0 : 0.0)
                             + 3.0 * fresh + selected
                             + 0.001 * label;
                choices.push_back({score, v});
            }
            int keep = min(branch, static_cast<int>(choices.size()));
            partial_sort(choices.begin(), choices.begin() + keep, choices.end(),
                         [](const auto& a, const auto& b) { return a.first > b.first; });
            for (int q = 0; q < keep; ++q) {
                BeamNode child = state;
                beam_add(p, child, choices[q].second);
                child.merit = beam_merit(p, child);
                next.push_back(std::move(child));
            }
        }
        int keep = min(width, static_cast<int>(next.size()));
        partial_sort(next.begin(), next.begin() + keep, next.end(),
                     [](const BeamNode& a, const BeamNode& b) { return a.merit > b.merit; });
        next.resize(keep);
        beam = std::move(next);
    }
    int best = 0;
    for (int i = 1; i < static_cast<int>(beam.size()); ++i)
        if (beam[i].cycles > beam[best].cycles) best = i;
    return std::move(beam[best].order);
}

class ExactSearch {
    struct Undo {
        bool merged = false;
        int parent_root = -1, child_root = -1;
        int old_size = 0, old_start = 0, old_finish = 0;
    };

    const Problem& p;
    Candidate& best;
    chrono::steady_clock::time_point deadline;
    vector<int> parent, sz, chain_start, chain_finish, order;
    vector<uint64_t> neighbor_mask;
    uint64_t nodes = 0;
    bool timed_out = false;

    int find_root(int x) const {
        while (parent[x] != x) x = parent[x];
        return x;
    }

    bool add_edge(int from, int to, Undo& undo) {
        int a = find_root(from), b = find_root(to);
        if (a == b) return true;
        int new_start = chain_start[a], new_finish = chain_finish[b];
        if (sz[a] < sz[b]) swap(a, b);
        undo = {true, a, b, sz[a], chain_start[a], chain_finish[a]};
        parent[b] = a;
        sz[a] += sz[b];
        chain_start[a] = new_start;
        chain_finish[a] = new_finish;
        return false;
    }

    void rollback(const Undo& u) {
        if (!u.merged) return;
        parent[u.child_root] = u.child_root;
        sz[u.parent_root] = u.old_size;
        chain_start[u.parent_root] = u.old_start;
        chain_finish[u.parent_root] = u.old_finish;
    }

    void dfs(uint64_t used, uint64_t frontier, int cycles) {
        if ((++nodes & 1023ULL) == 0 && chrono::steady_clock::now() >= deadline) {
            timed_out = true;
            return;
        }
        int k = static_cast<int>(order.size());
        if (k == p.n) {
            if (cycles > best.cycles) best = evaluate(p, order);
            return;
        }

        bool can_close_now = false;
        for (int v = 0; v < p.n; ++v) if (frontier & (1ULL << v)) {
            if (find_root(k) == find_root(p.rank_at[v])) {
                can_close_now = true;
                break;
            }
        }
        int upper = cycles + (p.n - k) - (can_close_now ? 0 : 1);
        if (upper <= best.cycles) return;

        vector<int> choices;
        for (int v = 0; v < p.n; ++v) if (frontier & (1ULL << v))
            choices.push_back(v);
        sort(choices.begin(), choices.end(), [&](int a, int b) {
            bool ca = find_root(k) == find_root(p.rank_at[a]);
            bool cb = find_root(k) == find_root(p.rank_at[b]);
            if (ca != cb) return ca > cb;
            return p.rank_at[a] < p.rank_at[b];
        });

        for (int v : choices) {
            Undo undo;
            bool closed = add_edge(k, p.rank_at[v], undo);
            order.push_back(v);
            uint64_t bit = 1ULL << v;
            uint64_t next_used = used | bit;
            uint64_t next_frontier = (frontier | neighbor_mask[v]) & ~next_used;
            dfs(next_used, next_frontier, cycles + closed);
            order.pop_back();
            rollback(undo);
            if (timed_out) return;
        }
    }

public:
    ExactSearch(const Problem& problem, Candidate& incumbent,
                chrono::steady_clock::time_point stop)
        : p(problem), best(incumbent), deadline(stop), parent(p.n), sz(p.n, 1),
          chain_start(p.n), chain_finish(p.n), neighbor_mask(p.n, 0) {
        iota(parent.begin(), parent.end(), 0);
        iota(chain_start.begin(), chain_start.end(), 0);
        iota(chain_finish.begin(), chain_finish.end(), 0);
        order.reserve(p.n);
        for (int v = 0; v < p.n; ++v)
            for (int u : p.nbr[v]) neighbor_mask[v] |= 1ULL << u;
    }

    bool run() {
        vector<int> roots(p.n);
        iota(roots.begin(), roots.end(), 0);
        stable_sort(roots.begin(), roots.end(), [&](int a, int b) {
            bool ar = a == p.pos[0], br = b == p.pos[0];
            if (ar != br) return ar;
            return a < b;
        });
        for (int root : roots) {
            Undo undo;
            bool closed = add_edge(0, p.rank_at[root], undo);
            order.push_back(root);
            dfs(1ULL << root, neighbor_mask[root] & ~(1ULL << root), closed);
            order.pop_back();
            rollback(undo);
            if (timed_out) break;
        }
        return !timed_out;
    }
};

static Candidate solve_for(const Problem& p, double seconds, SearchStats* stats = nullptr,
                           int iteration_limit = -1) {
    auto solve_start = chrono::steady_clock::now();
    auto elapsed_ms = [&]() {
        return 1000.0 * chrono::duration<double>(chrono::steady_clock::now() - solve_start).count();
    };
    Candidate best = evaluate(p, baseline_order(p));
    adjacent_improve(p, best);
    if (stats) stats->fallback_swaps = p.n - best.cycles;
    if (p.n <= 1 || seconds <= 0.0) return best;

    mt19937_64 rng(matrix_seed(p));
    const vector<RepairWeights> repair_deterministic = {
        {1000, 10, 1, 0, 0, 0, 16},
        {100, 5, 3, -1, 2, 0, 24},
        {30, 15, -1, 2, -2, 0, 10},
        {3000, 0, 0, 0, 0, 0, 8},
        {4.092705, 44.818952, 6.683821, -0.545920, 4.979314, 0, 64},
        {4.142114, 78.174754, 9.560245, 0.256275, 7.571936, 0, 48},
        {5.527978, 79.822470, 5.858314, -5.044624, 5.824673, 0, 12},
    };
    for (int repair_index = 0;
         repair_index < static_cast<int>(repair_deterministic.size()); ++repair_index) {
        // The learned specialists were trained on large, high-frontier boards;
        // below 60 cells they add cost and can perturb a stronger small-board
        // search trajectory without a consistent cycle benefit.
        if (repair_index >= 4 && p.n < 60) continue;
        const RepairWeights& w = repair_deterministic[repair_index];
        Candidate c = evaluate(p, repair_order(p, rng, w));
        adjacent_improve(p, c);
        if (better(c, best)) {
            best = std::move(c);
            if (stats) ++stats->deterministic_wins;
        }
    }
    const vector<Weights> deterministic = {
        {10000, 1000, 80, 4, 1, 0.5, 0, 0, 0, 16},
        {10000, 300, 30, 2, 2, 0.2, 2, -1, 0, 10},
        {10000, 100, 10, 1, -1, 1.5, -2, 2, 0, 24},
        {10000, 3000, 150, 0, 0.5, 0, 1, -2, 0, 8},
        {300, 100, 20, 3, 1, 0.5, 0, 0, 0, 16},
    };
    for (const Weights& w : deterministic) {
        Candidate c = evaluate(p, greedy_order(p, rng, w, p.pos[0]));
        adjacent_improve(p, c);
        if (better(c, best)) {
            best = std::move(c);
            if (stats) ++stats->deterministic_wins;
        }
    }
    if (stats) {
        stats->deterministic_swaps = p.n - best.cycles;
        stats->deterministic_ms = elapsed_ms();
    }

    auto raw_deadline = chrono::steady_clock::now() + chrono::duration<double>(seconds);
    auto deadline = chrono::time_point_cast<chrono::steady_clock::duration>(raw_deadline);
    int cycles_before_local = best.cycles;
    transposition_improve(p, best, rng, 30 * p.n, deadline);
    if (stats) {
        stats->local_wins += best.cycles - cycles_before_local;
        stats->post_local_swaps = p.n - best.cycles;
        stats->local_ms = elapsed_ms();
    }
    if (chrono::steady_clock::now() < deadline) {
        auto beam_raw_stop = chrono::steady_clock::now()
                           + chrono::duration<double>(min(2.0, seconds * 0.12));
        auto beam_stop = min(deadline,
            chrono::time_point_cast<chrono::steady_clock::duration>(beam_raw_stop));
        auto consider_beam = [&](vector<int> order) {
            if (order.empty()) return;
            Candidate c = evaluate(p, std::move(order));
            adjacent_improve(p, c);
            if (better(c, best)) {
                best = std::move(c);
                if (stats) {
                    ++stats->beam_wins;
                    stats->last_improvement_ms = elapsed_ms();
                }
            }
        };
        if (p.n > 225) {
            int width = seconds >= 2.0 ? 64 : 32;
            consider_beam(repair_beam_order(
                p, repair_deterministic[4], width, 3, beam_stop));
            if (seconds >= 2.0 && chrono::steady_clock::now() < beam_stop)
                consider_beam(repair_beam_order(
                    p, repair_deterministic[5], 32, 4, beam_stop));
            if (seconds >= 2.0 && chrono::steady_clock::now() < beam_stop)
                consider_beam(repair_beam_order(
                    p, repair_deterministic[6], 32, 3, beam_stop));
        } else {
            consider_beam(beam_order(p, 64, 4, beam_stop));
        }
    }
    if (stats) {
        stats->post_beam_swaps = p.n - best.cycles;
        stats->beam_ms = elapsed_ms();
    }
    if (p.n <= 36 && chrono::steady_clock::now() < deadline) {
        auto exact_raw_stop = chrono::steady_clock::now()
                            + chrono::duration<double>(min(5.0, seconds * 0.20));
        auto exact_stop = min(deadline,
            chrono::time_point_cast<chrono::steady_clock::duration>(exact_raw_stop));
        ExactSearch exact(p, best, exact_stop);
        // Even after proving the primary optimum, keep the remaining anytime
        // budget: other optimal-swap targets can improve maxWorstDistance.
        int cycles_before_exact = best.cycles;
        exact.run();
        if (stats && best.cycles > cycles_before_exact) {
            stats->exact_wins += best.cycles - cycles_before_exact;
            stats->last_improvement_ms = elapsed_ms();
        }
    }
    if (stats) {
        stats->post_exact_swaps = p.n - best.cycles;
        stats->exact_ms = elapsed_ms();
    }
    int iteration = 0;
    while (chrono::steady_clock::now() < deadline &&
           (iteration_limit < 0 || iteration < iteration_limit)) {
        if (p.n > 225 && (iteration & 3) == 0) {
            RepairWeights rw = repair_deterministic[4 + rng() % 3];
            rw.expose_displaced *= 0.3 + double(rng() % 2700) / 1000.0;
            rw.expose_upcoming *= 0.3 + double(rng() % 2400) / 1000.0;
            rw.new_frontier += double(int(rng() % 601) - 300) / 100.0;
            rw.selected_neighbors += double(int(rng() % 401) - 200) / 100.0;
            rw.delay += double(int(rng() % 601) - 300) / 100.0;
            int width = (iteration & 4) ? 16 : 8;
            vector<int> order = repair_beam_order(p, rw, width, 3, deadline);
            if (!order.empty()) {
                Candidate c = evaluate(p, std::move(order));
                adjacent_improve(p, c);
                if (better(c, best)) {
                    best = std::move(c);
                    if (stats) {
                        ++stats->random_beam_wins;
                        stats->last_improvement_ms = elapsed_ms();
                    }
                    int before = best.cycles;
                    transposition_improve(p, best, rng, 12 * p.n, deadline);
                    if (stats) stats->local_wins += best.cycles - before;
                }
            }
            ++iteration;
            continue;
        }
        if ((iteration & 1) == 1) {
            int repair_family = p.n >= 60 ? int(rng() % 16) : -1;
            RepairWeights rw = p.n >= 60
                ? repair_deterministic[repair_family < 14 ? 6 : 4 + (repair_family & 1)]
                : repair_deterministic[rng() % 4];
            rw.noise = 2.0 + double(rng() % 10000) / 100.0;
            rw.expose_displaced *= 0.2 + double(rng() % 3000) / 1000.0;
            rw.expose_upcoming *= double(rng() % 3000) / 1000.0;
            rw.new_frontier += double(int(rng() % 801) - 400) / 100.0;
            rw.delay += double(int(rng() % 801) - 400) / 50.0;
            Candidate c = evaluate(p, repair_order(p, rng, rw));
            adjacent_improve(p, c);
            if (better(c, best)) {
                best = std::move(c);
                if (stats) {
                    ++stats->random_repair_wins;
                    stats->last_improvement_ms = elapsed_ms();
                }
                int before = best.cycles;
                transposition_improve(p, best, rng, 12 * p.n, deadline);
                if (stats) stats->local_wins += best.cycles - before;
            }
            ++iteration;
            continue;
        }
        Weights w = deterministic[rng() % deterministic.size()];
        w.noise = 5.0 + double(rng() % 4000) / 100.0;
        w.reserved *= 0.5 + double(rng() % 2000) / 1000.0;
        w.ready *= 0.5 + double(rng() % 2000) / 1000.0;
        w.expose *= double(rng() % 2000) / 1000.0;
        w.new_frontier += double(int(rng() % 601) - 300) / 100.0;
        w.future_rank += double(int(rng() % 601) - 300) / 50.0;

        vector<int> order;
        bool used_destroy_repair = (iteration & 3) != 0 && p.n >= 8;
        if (used_destroy_repair) {
            int hi = max(2, p.n - 2);
            int cut = 1 + int(rng() % hi);
            // Bias toward large retained prefixes, with occasional broad restart.
            if (rng() & 1) cut = max(1, p.n - 1 - int(rng() % max(2, p.n / 3)));
            order = greedy_order(p, rng, w, best.order[0], &best.order, cut);
        } else {
            ++iteration;
            int root = (iteration % 17 == 0) ? int(rng() % p.n) : p.pos[0];
            order = greedy_order(p, rng, w, root);
        }
        if (used_destroy_repair) ++iteration;
        Candidate c = evaluate(p, std::move(order));
        adjacent_improve(p, c);
        if (better(c, best)) {
            best = std::move(c);
            if (stats) {
                if (used_destroy_repair) ++stats->destroy_repair_wins;
                else ++stats->random_restart_wins;
                stats->last_improvement_ms = elapsed_ms();
            }
            int before = best.cycles;
            transposition_improve(p, best, rng, 12 * p.n, deadline);
            if (stats) stats->local_wins += best.cycles - before;
        }
    }
    if (stats) {
        stats->iterations = iteration;
        stats->final_swaps = p.n - best.cycles;
    }
    return best;
}

[[maybe_unused]] static bool connected_order(const Problem& p, const vector<int>& order) {
    if (static_cast<int>(order.size()) != p.n) return false;
    vector<unsigned char> used(p.n, 0);
    for (int k = 0; k < p.n; ++k) {
        int v = order[k];
        if (v < 0 || v >= p.n || used[v]) return false;
        if (k > 0) {
            bool touches = false;
            for (int u : p.nbr[v]) touches |= used[u];
            if (!touches) return false;
        }
        used[v] = 1;
    }
    return true;
}

static vector<vector<int>> swaps_for_order(const Problem& p,
                                            const vector<int>& order) {
    vector<int> desired(p.n), cur = p.rank_at, location(p.n);
    for (int k = 0; k < p.n; ++k) desired[order[k]] = k;
    for (int v = 0; v < p.n; ++v) location[cur[v]] = v;

    vector<vector<int>> swaps;
    swaps.reserve(p.n);
    for (int v = 0; v < p.n; ++v) {
        if (cur[v] == desired[v]) continue;
        int u = location[desired[v]];
        int rv = cur[v], ru = cur[u];
        swap(cur[v], cur[u]);
        location[rv] = u;
        location[ru] = v;
        swaps.push_back({v / p.cols, v % p.cols, u / p.cols, u % p.cols});
    }
    return swaps;
}

}  // namespace hollan

vector<vector<int>> get_swaps(vector<vector<int>> matrix) {
    hollan::Problem problem(matrix);
    // The official allowance is 120 wall-clock seconds. Ninety-eight seconds
    // leaves ample time for wrapper I/O, scheduling jitter, and submission.
    hollan::Candidate best = hollan::solve_for(problem, 98.0);
    return hollan::swaps_for_order(problem, best.order);
}
