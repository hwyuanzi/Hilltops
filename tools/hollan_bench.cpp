#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "../sample_bots/hollan_bot.cpp"

using namespace std;

static pair<bool, int> authoritative_eval(const vector<vector<int>>& a) {
    int R = static_cast<int>(a.size()), C = static_cast<int>(a[0].size());
    int n = R * C;
    vector<int> order(n), reachable(n, 0), worst(n, -1);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int x, int y) {
        return a[x / C][x % C] < a[y / C][y % C];
    });
    reachable[order[0]] = 1;
    worst[order[0]] = 0;
    int answer = 0;
    for (int k = 1; k < n; ++k) {
        int v = order[k], r = v / C, c = v % C;
        static const int dr[4] = {-1, 1, 0, 0};
        static const int dc[4] = {0, 0, -1, 1};
        for (int d = 0; d < 4; ++d) {
            int nr = r + dr[d], nc = c + dc[d];
            if (0 <= nr && nr < R && 0 <= nc && nc < C) {
                int u = nr * C + nc;
                if (reachable[u] && a[nr][nc] < a[r][c]) {
                    reachable[v] = 1;
                    worst[v] = max(worst[v], worst[u] + 1);
                }
            }
        }
        if (!reachable[v]) return {false, -1};
        answer = max(answer, worst[v]);
    }
    return {true, answer};
}

static vector<vector<int>> apply_swaps(vector<vector<int>> a,
                                       const vector<vector<int>>& swaps) {
    for (const auto& s : swaps) {
        assert(s.size() == 4);
        assert(0 <= s[0] && s[0] < static_cast<int>(a.size()));
        assert(0 <= s[2] && s[2] < static_cast<int>(a.size()));
        assert(0 <= s[1] && s[1] < static_cast<int>(a[0].size()));
        assert(0 <= s[3] && s[3] < static_cast<int>(a[0].size()));
        swap(a[s[0]][s[1]], a[s[2]][s[3]]);
    }
    return a;
}

int main(int argc, char** argv) {
    int trials = argc > 1 ? stoi(argv[1]) : 200;
    mt19937_64 rng(0x48494c4c544f5053ULL);
    const vector<pair<int, int>> sizes = {
        {3, 3}, {4, 4}, {5, 5}, {8, 8}, {10, 10}, {15, 15}, {20, 20}};
    for (auto [R, C] : sizes) {
        long long total_swaps = 0, total_worst = 0;
        int min_swaps = R * C, max_swaps = 0;
        auto start = chrono::steady_clock::now();
        for (int t = 0; t < trials; ++t) {
            vector<int> values(R * C);
            iota(values.begin(), values.end(), 1);
            shuffle(values.begin(), values.end(), rng);
            vector<vector<int>> matrix(R, vector<int>(C));
            for (int i = 0; i < R * C; ++i) matrix[i / C][i % C] = values[i];
            auto swaps = get_swaps(matrix);
            auto final_matrix = apply_swaps(matrix, swaps);
            auto [ok, worst] = authoritative_eval(final_matrix);
            if (!ok) {
                cerr << "INVALID " << R << 'x' << C << " trial " << t << '\n';
                return 2;
            }
            hollan::Problem p(matrix);
            auto target = hollan::baseline_order(p);
            assert(static_cast<int>(swaps.size()) == p.n - hollan::cycle_count(p, target));
            total_swaps += swaps.size();
            total_worst += worst;
            min_swaps = min(min_swaps, static_cast<int>(swaps.size()));
            max_swaps = max(max_swaps, static_cast<int>(swaps.size()));
        }
        double seconds = chrono::duration<double>(chrono::steady_clock::now() - start).count();
        cout << R << 'x' << C
             << " valid=" << trials << '/' << trials
             << " swaps_avg=" << double(total_swaps) / trials
             << " swaps_range=" << min_swaps << '-' << max_swaps
             << " worst_avg=" << double(total_worst) / trials
             << " ms_per=" << 1000.0 * seconds / trials << '\n';
    }
}
