#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "../sample_bots/hollan_bot.cpp"

using namespace std;

struct Case {
    uint64_t seed;
    int rows, cols;
    string family;
    vector<vector<int>> matrix;
};

static uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static vector<vector<int>> random_matrix(int rows, int cols, uint64_t seed) {
    vector<int> a(rows * cols);
    iota(a.begin(), a.end(), 1);
    uint64_t state = seed;
    for (int i = static_cast<int>(a.size()) - 1; i > 0; --i)
        swap(a[i], a[splitmix64(state) % static_cast<uint64_t>(i + 1)]);
    vector<vector<int>> matrix(rows, vector<int>(cols));
    for (int i = 0; i < rows * cols; ++i) matrix[i / cols][i % cols] = a[i];
    return matrix;
}

static vector<vector<int>> pattern_matrix(int rows, int cols, int kind) {
    int n = rows * cols;
    vector<int> placement(n);
    if (kind == 0) {
        iota(placement.begin(), placement.end(), 1);
    } else if (kind == 1) {
        for (int i = 0; i < n; ++i) placement[i] = n - i;
    } else if (kind == 2) {
        int rank = 1;
        for (int r = 0; r < rows; ++r) {
            if ((r & 1) == 0) for (int c = 0; c < cols; ++c) placement[r * cols + c] = rank++;
            else for (int c = cols - 1; c >= 0; --c) placement[r * cols + c] = rank++;
        }
    } else if (kind == 3) {
        vector<int> cells(n);
        iota(cells.begin(), cells.end(), 0);
        sort(cells.begin(), cells.end(), [&](int a, int b) {
            int ar = a / cols, ac = a % cols, br = b / cols, bc = b % cols;
            int da = min(min(ar, rows - 1 - ar), min(ac, cols - 1 - ac));
            int db = min(min(br, rows - 1 - br), min(bc, cols - 1 - bc));
            if (da != db) return da < db;
            return a < b;
        });
        for (int k = 0; k < n; ++k) placement[cells[k]] = k + 1;
    } else {
        vector<int> low, high;
        int lo = 1, hi = n;
        for (int i = 0; i < n; ++i) placement[i] = ((i / cols + i % cols) & 1) ? hi-- : lo++;
    }
    vector<vector<int>> matrix(rows, vector<int>(cols));
    for (int i = 0; i < n; ++i) matrix[i / cols][i % cols] = placement[i];
    return matrix;
}

static void add_random(vector<Case>& out, int rows, int cols, int count,
                       uint64_t base, const string& family) {
    for (int i = 0; i < count; ++i) {
        uint64_t seed = base + static_cast<uint64_t>(i);
        out.push_back({seed, rows, cols, family, random_matrix(rows, cols, seed)});
    }
}

static vector<Case> make_corpus(const string& name) {
    vector<Case> out;
    uint64_t base = name == "heldout" ? 0xd1b54a32d192ed03ULL : 0x243f6a8885a308d3ULL;
    if (name == "full" || name == "train") {
        add_random(out, 5, 5, 1000, base + 0x10000, "random");
        add_random(out, 8, 8, 500, base + 0x20000, "random");
        add_random(out, 10, 10, 500, base + 0x30000, "random");
        add_random(out, 15, 15, 200, base + 0x40000, "random");
        add_random(out, 20, 20, 200, base + 0x50000, "random");
        add_random(out, 1, 20, 100, base + 0x60000, "thin");
        add_random(out, 2, 20, 100, base + 0x70000, "thin");
        add_random(out, 3, 20, 100, base + 0x80000, "thin");
        add_random(out, 4, 10, 100, base + 0x90000, "rectangle");
        add_random(out, 6, 12, 100, base + 0xa0000, "rectangle");
        add_random(out, 8, 15, 100, base + 0xb0000, "rectangle");
        for (auto [r, c] : vector<pair<int, int>>{{5, 5}, {8, 8}, {10, 10}, {15, 15}, {20, 20}})
            for (int kind = 0; kind < 5; ++kind)
                out.push_back({0xf0000000ULL + static_cast<uint64_t>(1000 * r + kind), r, c,
                               "pattern" + to_string(kind), pattern_matrix(r, c, kind)});
    } else if (name == "heldout") {
        add_random(out, 5, 5, 200, base + 0x10000, "random");
        add_random(out, 8, 8, 100, base + 0x20000, "random");
        add_random(out, 10, 10, 100, base + 0x30000, "random");
        add_random(out, 15, 15, 50, base + 0x40000, "random");
        add_random(out, 20, 20, 50, base + 0x50000, "random");
        add_random(out, 2, 20, 30, base + 0x60000, "thin");
        add_random(out, 3, 20, 30, base + 0x70000, "thin");
        add_random(out, 7, 13, 30, base + 0x80000, "rectangle");
        add_random(out, 11, 17, 30, base + 0x90000, "rectangle");
    } else if (name == "20x20") {
        add_random(out, 20, 20, 50, base + 0x50000, "random");
    } else {
        cerr << "unknown corpus: " << name << '\n';
        exit(2);
    }
    return out;
}

static pair<bool, int> independent_eval(const vector<vector<int>>& a) {
    int rows = static_cast<int>(a.size()), cols = static_cast<int>(a[0].size());
    int n = rows * cols;
    vector<int> order(n), dist(n, -1);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int x, int y) {
        int vx = a[x / cols][x % cols], vy = a[y / cols][y % cols];
        return vx != vy ? vx < vy : x < y;
    });
    dist[order[0]] = 0;
    int worst = 0;
    for (int k = 1; k < n; ++k) {
        int v = order[k], r = v / cols, c = v % cols;
        static const int dr[] = {-1, 1, 0, 0};
        static const int dc[] = {0, 0, -1, 1};
        for (int d = 0; d < 4; ++d) {
            int rr = r + dr[d], cc = c + dc[d];
            if (0 <= rr && rr < rows && 0 <= cc && cc < cols) {
                int u = rr * cols + cc;
                if (a[rr][cc] < a[r][c] && dist[u] >= 0)
                    dist[v] = max(dist[v], dist[u] + 1);
            }
        }
        if (dist[v] < 0) return {false, -1};
        worst = max(worst, dist[v]);
    }
    return {true, worst};
}

static int independent_cycles(const hollan::Problem& p, const vector<int>& order) {
    vector<int> permutation(p.n), seen(p.n, 0);
    for (int k = 0; k < p.n; ++k) permutation[k] = p.rank_at[order[k]];
    int cycles = 0;
    for (int s = 0; s < p.n; ++s) if (!seen[s]) {
        ++cycles;
        for (int v = s; !seen[v]; v = permutation[v]) seen[v] = 1;
    }
    return cycles;
}

int main(int argc, char** argv) {
    string corpus_name = argc > 1 ? argv[1] : "heldout";
    double seconds = argc > 2 ? stod(argv[2]) : 0.02;
    string output = argc > 3 ? argv[3] : "hollan_results.csv";
    int limit = argc > 4 ? stoi(argv[4]) : -1;
    vector<Case> corpus = make_corpus(corpus_name);
    if (limit >= 0 && limit < static_cast<int>(corpus.size())) corpus.resize(limit);
    ofstream csv(output);
    if (!csv) { cerr << "cannot write " << output << '\n'; return 2; }
    csv << "index,seed,rows,cols,family,seconds,swaps,worst,runtime_ms,valid\n";
    long long total_swaps = 0, total_worst = 0;
    double total_ms = 0.0, max_ms = 0.0;
    int valid_count = 0;
    for (int index = 0; index < static_cast<int>(corpus.size()); ++index) {
        const Case& tc = corpus[index];
        hollan::Problem p(tc.matrix);
        auto start = chrono::steady_clock::now();
        hollan::Candidate candidate = hollan::solve_for(p, seconds);
        vector<vector<int>> swaps = hollan::swaps_for_order(p, candidate.order);
        double elapsed_ms = 1000.0 * chrono::duration<double>(chrono::steady_clock::now() - start).count();

        vector<vector<int>> final_matrix = tc.matrix;
        bool legal = static_cast<int>(candidate.order.size()) == p.n;
        for (const vector<int>& s : swaps) {
            if (s.size() != 4 || s[0] < 0 || s[0] >= tc.rows || s[2] < 0 || s[2] >= tc.rows ||
                s[1] < 0 || s[1] >= tc.cols || s[3] < 0 || s[3] >= tc.cols) {
                legal = false;
                break;
            }
            swap(final_matrix[s[0]][s[1]], final_matrix[s[2]][s[3]]);
        }
        vector<vector<int>> intended(tc.rows, vector<int>(tc.cols));
        if (legal) for (int k = 0; k < p.n; ++k)
            intended[candidate.order[k] / tc.cols][candidate.order[k] % tc.cols] = k + 1;
        bool exact_target = legal && final_matrix == intended;
        auto [hilltop, worst] = independent_eval(final_matrix);
        int cycles = legal ? independent_cycles(p, candidate.order) : -1;
        bool valid = legal && exact_target && hilltop &&
                     static_cast<int>(swaps.size()) == p.n - cycles &&
                     candidate.cycles == cycles && candidate.worst == worst;
        valid_count += valid;
        total_swaps += swaps.size();
        total_worst += worst;
        total_ms += elapsed_ms;
        max_ms = max(max_ms, elapsed_ms);
        csv << index << ',' << tc.seed << ',' << tc.rows << ',' << tc.cols << ',' << tc.family << ','
            << setprecision(12) << seconds << ',' << swaps.size() << ',' << worst << ','
            << elapsed_ms << ',' << (valid ? 1 : 0) << '\n';
        if (!valid) {
            cerr << "INVALID index=" << index << " seed=" << tc.seed << " size="
                 << tc.rows << 'x' << tc.cols << '\n';
            return 3;
        }
        if ((index + 1) % 100 == 0)
            cerr << "completed " << index + 1 << '/' << corpus.size() << '\n';
    }
    cout << "corpus=" << corpus_name << " cases=" << corpus.size()
         << " valid=" << valid_count << '/' << corpus.size()
         << " swaps_avg=" << fixed << setprecision(4) << double(total_swaps) / corpus.size()
         << " worst_avg=" << double(total_worst) / corpus.size()
         << " runtime_ms_avg=" << total_ms / corpus.size()
         << " runtime_ms_max=" << max_ms << " output=" << output << '\n';
}
