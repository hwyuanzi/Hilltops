#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "../sample_bots/hollan_bot.cpp"

using namespace std;

static uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static vector<vector<int>> matrix_for(uint64_t seed) {
    vector<int> a(400);
    iota(a.begin(), a.end(), 1);
    uint64_t state = seed;
    for (int i = 399; i > 0; --i) swap(a[i], a[splitmix64(state) % uint64_t(i + 1)]);
    vector<vector<int>> matrix(20, vector<int>(20));
    for (int i = 0; i < 400; ++i) matrix[i / 20][i % 20] = a[i];
    return matrix;
}

struct Item {
    hollan::RepairWeights w;
    int saved = 0, wins = 0, max_saved = 0;
};

static double uniform(mt19937_64& rng, double lo, double hi) {
    return lo + (hi - lo) * double(rng() >> 11) * (1.0 / 9007199254740992.0);
}

int main(int argc, char** argv) {
    int train_count = argc > 1 ? stoi(argv[1]) : 100;
    int parameter_count = argc > 2 ? stoi(argv[2]) : 1000;
    uint64_t base_seed = 0x243f6a8885a308d3ULL + 0x50000ULL;
    vector<hollan::Problem> problems;
    problems.reserve(train_count);
    for (int i = 0; i < train_count; ++i) problems.emplace_back(matrix_for(base_seed + i));

    const vector<hollan::RepairWeights> repairs = {
        {1000, 10, 1, 0, 0, 0, 16},
        {100, 5, 3, -1, 2, 0, 24},
        {30, 15, -1, 2, -2, 0, 10},
        {3000, 0, 0, 0, 0, 0, 8},
    };
    const vector<hollan::Weights> greedies = {
        {10000, 1000, 80, 4, 1, 0.5, 0, 0, 0, 16},
        {10000, 300, 30, 2, 2, 0.2, 2, -1, 0, 10},
        {10000, 100, 10, 1, -1, 1.5, -2, 2, 0, 24},
        {10000, 3000, 150, 0, 0.5, 0, 1, -2, 0, 8},
        {300, 100, 20, 3, 1, 0.5, 0, 0, 0, 16},
    };
    vector<int> base_cycles(train_count, 0);
    for (int i = 0; i < train_count; ++i) {
        const auto& p = problems[i];
        mt19937_64 rng(hollan::matrix_seed(p));
        auto consider = [&](vector<int> order) {
            hollan::Candidate c = hollan::evaluate(p, std::move(order));
            hollan::adjacent_improve(p, c);
            base_cycles[i] = max(base_cycles[i], c.cycles);
        };
        consider(hollan::baseline_order(p));
        for (const auto& w : repairs) consider(hollan::repair_order(p, rng, w));
        for (const auto& w : greedies) consider(hollan::greedy_order(p, rng, w, p.pos[0]));
    }

    mt19937_64 search_rng(0x504152414d534541ULL);
    vector<Item> items;
    items.reserve(parameter_count);
    const vector<int> looks = {4, 8, 12, 16, 24, 32, 48, 64};
    for (int q = 0; q < parameter_count; ++q) {
        hollan::RepairWeights w;
        w.expose_displaced = exp(uniform(search_rng, log(3.0), log(10000.0)));
        w.expose_upcoming = uniform(search_rng, -30.0, 100.0);
        w.new_frontier = uniform(search_rng, -12.0, 12.0);
        w.selected_neighbors = uniform(search_rng, -8.0, 8.0);
        w.delay = uniform(search_rng, -30.0, 30.0);
        w.noise = 0.0;
        w.lookahead = looks[search_rng() % looks.size()];
        Item item{w};
        for (int i = 0; i < train_count; ++i) {
            const auto& p = problems[i];
            mt19937_64 rng(hollan::matrix_seed(p));
            hollan::Candidate c = hollan::evaluate(p, hollan::repair_order(p, rng, w));
            hollan::adjacent_improve(p, c);
            int saved = max(0, c.cycles - base_cycles[i]);
            item.saved += saved;
            item.wins += saved > 0;
            item.max_saved = max(item.max_saved, saved);
        }
        items.push_back(item);
        if ((q + 1) % 100 == 0) cerr << "parameters " << q + 1 << '/' << parameter_count << '\n';
    }
    sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        if (a.saved != b.saved) return a.saved > b.saved;
        return a.wins > b.wins;
    });
    cout << fixed << setprecision(6);
    for (int i = 0; i < min(20, static_cast<int>(items.size())); ++i) {
        const Item& x = items[i];
        cout << i << " saved=" << x.saved << " wins=" << x.wins << " max=" << x.max_saved
             << " weights={" << x.w.expose_displaced << ',' << x.w.expose_upcoming << ','
             << x.w.new_frontier << ',' << x.w.selected_neighbors << ',' << x.w.delay << ",0,"
             << x.w.lookahead << "}\n";
    }
}
