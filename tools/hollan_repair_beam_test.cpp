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
    int trials = argc > 1 ? stoi(argv[1]) : 50;
    int weight_index = argc > 2 ? stoi(argv[2]) : 0;
    const vector<hollan::RepairWeights> weights = {
        {4.092705, 44.818952, 6.683821, -0.545920, 4.979314, 0, 64},
        {4.142114, 78.174754, 9.560245, 0.256275, 7.571936, 0, 48},
        {5.527978, 79.822470, 5.858314, -5.044624, 5.824673, 0, 12},
    };
    const hollan::RepairWeights& w = weights.at(weight_index);
    const vector<pair<int, int>> configs = {{8, 3}, {16, 3}, {32, 3}, {32, 4}, {64, 3}};
    vector<long long> swaps(configs.size()), wins(configs.size());
    vector<double> elapsed(configs.size());
    long long greedy_swaps = 0, combo_swaps = 0, old_combo_swaps = 0;
    double combo_elapsed = 0.0;
    for (int t = 0; t < trials; ++t) {
        hollan::Problem p(matrix_for(0x243f6a8885a308d3ULL + 0x50000ULL + t));
        mt19937_64 rng(hollan::matrix_seed(p));
        hollan::Candidate greedy = hollan::evaluate(p, hollan::repair_order(p, rng, w));
        hollan::adjacent_improve(p, greedy);
        greedy_swaps += p.n - greedy.cycles;
        for (int q = 0; q < static_cast<int>(configs.size()); ++q) {
            auto start = chrono::steady_clock::now();
            auto deadline = start + chrono::seconds(30);
            vector<int> order = hollan::repair_beam_order(
                p, w, configs[q].first, configs[q].second, deadline);
            elapsed[q] += chrono::duration<double, milli>(chrono::steady_clock::now() - start).count();
            if (order.empty() || !hollan::connected_order(p, order)) return 2;
            hollan::Candidate c = hollan::evaluate(p, std::move(order));
            hollan::adjacent_improve(p, c);
            swaps[q] += p.n - c.cycles;
            wins[q] += c.cycles > greedy.cycles;
        }
        auto combo_start = chrono::steady_clock::now();
        int combo_best = 0;
        const vector<pair<int, int>> combo = {{64, 3}, {32, 4}, {32, 3}};
        for (int wi = 0; wi < 3; ++wi) {
            auto deadline = chrono::steady_clock::now() + chrono::seconds(30);
            vector<int> order = hollan::repair_beam_order(
                p, weights[wi], combo[wi].first, combo[wi].second, deadline);
            hollan::Candidate c = hollan::evaluate(p, std::move(order));
            hollan::adjacent_improve(p, c);
            combo_best = max(combo_best, c.cycles);
        }
        combo_swaps += p.n - combo_best;
        int old_combo_best = 0;
        for (int wi = 0; wi < 3; ++wi) {
            auto deadline = chrono::steady_clock::now() + chrono::seconds(30);
            vector<int> order = hollan::repair_beam_order(
                p, weights[wi], combo[wi].first, combo[wi].second, deadline, false);
            hollan::Candidate c = hollan::evaluate(p, std::move(order));
            hollan::adjacent_improve(p, c);
            old_combo_best = max(old_combo_best, c.cycles);
        }
        old_combo_swaps += p.n - old_combo_best;
        combo_elapsed += chrono::duration<double, milli>(
            chrono::steady_clock::now() - combo_start).count();
        cerr << "beam test " << t + 1 << '/' << trials << '\n';
    }
    cout << "weight=" << weight_index << " greedy avg=" << double(greedy_swaps) / trials << '\n';
    for (int q = 0; q < static_cast<int>(configs.size()); ++q)
        cout << "width=" << configs[q].first << " branch=" << configs[q].second
             << " avg=" << double(swaps[q]) / trials << " wins=" << wins[q]
             << " ms=" << elapsed[q] / trials << '\n';
    cout << "combo avg=" << double(combo_swaps) / trials
         << " old_combo=" << double(old_combo_swaps) / trials
         << " ms=" << combo_elapsed / trials << '\n';
}
