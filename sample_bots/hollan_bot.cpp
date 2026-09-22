#include <algorithm>
#include <chrono>
#include <cstdint>
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

static bool connected_order(const Problem& p, const vector<int>& order) {
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
    vector<int> order = hollan::baseline_order(problem);
    return hollan::swaps_for_order(problem, order);
}
