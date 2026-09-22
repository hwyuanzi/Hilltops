#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
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

int main(int argc, char** argv) {
    int trials = argc > 1 ? stoi(argv[1]) : 30;
    int restarts = argc > 2 ? stoi(argv[2]) : 50;
    int width = argc > 3 ? stoi(argv[3]) : 8;
    const vector<hollan::RepairWeights> bases = {
        {4.092705, 44.818952, 6.683821, -0.545920, 4.979314, 0, 64},
        {4.142114, 78.174754, 9.560245, 0.256275, 7.571936, 0, 48},
        {5.527978, 79.822470, 5.858314, -5.044624, 5.824673, 0, 12},
    };
    long long before = 0, after = 0;
    int wins = 0, max_saved = 0;
    double search_ms = 0.0;
    for (int t = 0; t < trials; ++t) {
        hollan::Problem p(matrix_for(0x243f6a8885a308d3ULL + 0x50000ULL + t));
        hollan::Candidate best = hollan::solve_for(p, 2.0, nullptr, 0);
        int old = p.n - best.cycles;
        mt19937_64 rng(hollan::matrix_seed(p) ^ 0x52414e444245414dULL);
        auto start = chrono::steady_clock::now();
        for (int q = 0; q < restarts; ++q) {
            hollan::RepairWeights w = bases[rng() % bases.size()];
            w.expose_displaced *= 0.3 + double(rng() % 2700) / 1000.0;
            w.expose_upcoming *= 0.3 + double(rng() % 2400) / 1000.0;
            w.new_frontier += double(int(rng() % 601) - 300) / 100.0;
            w.selected_neighbors += double(int(rng() % 401) - 200) / 100.0;
            w.delay += double(int(rng() % 601) - 300) / 100.0;
            auto deadline = chrono::steady_clock::now() + chrono::seconds(30);
            vector<int> order = hollan::repair_beam_order(p, w, width, 3, deadline);
            hollan::Candidate c = hollan::evaluate(p, std::move(order));
            hollan::adjacent_improve(p, c);
            if (hollan::better(c, best)) best = std::move(c);
        }
        search_ms += chrono::duration<double, milli>(chrono::steady_clock::now() - start).count();
        int now = p.n - best.cycles;
        before += old;
        after += now;
        wins += now < old;
        max_saved = max(max_saved, old - now);
        cerr << "random beam " << t + 1 << '/' << trials << ' ' << old << " -> " << now << '\n';
    }
    cout << "width=" << width << " restarts=" << restarts
         << " before=" << double(before) / trials << " after=" << double(after) / trials
         << " wins=" << wins << '/' << trials << " max_saved=" << max_saved
         << " ms=" << search_ms / trials << '\n';
}
