#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "../sample_bots/hollan_bot.cpp"

using namespace std;

int main(int argc, char** argv) {
    int trials = argc > 1 ? stoi(argv[1]) : 100;
    int side = argc > 2 ? stoi(argv[2]) : 20;
    mt19937_64 data_rng(0x48494c4c544f5053ULL);
    vector<pair<string, long long>> totals;
    auto add = [&](int index, const string& name, int swaps) {
        if (index >= static_cast<int>(totals.size())) totals.push_back({name, 0});
        totals[index].second += swaps;
    };
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

    for (int t = 0; t < trials; ++t) {
        vector<int> values(side * side);
        iota(values.begin(), values.end(), 1);
        shuffle(values.begin(), values.end(), data_rng);
        vector<vector<int>> matrix(side, vector<int>(side));
        for (int i = 0; i < side * side; ++i) matrix[i / side][i % side] = values[i];
        hollan::Problem p(matrix);
        mt19937_64 rng(hollan::matrix_seed(p));
        int index = 0;
        auto score = [&](vector<int> order, const string& name) {
            hollan::Candidate c = hollan::evaluate(p, std::move(order));
            hollan::adjacent_improve(p, c);
            add(index++, name, p.n - c.cycles);
        };
        score(hollan::baseline_order(p), "baseline");
        for (int i = 0; i < static_cast<int>(repairs.size()); ++i)
            score(hollan::repair_order(p, rng, repairs[i]), "repair" + to_string(i));
        for (int i = 0; i < static_cast<int>(greedies.size()); ++i)
            score(hollan::greedy_order(p, rng, greedies[i], p.pos[0]), "greedy" + to_string(i));
        auto stop = chrono::steady_clock::now() + chrono::seconds(60);
        score(hollan::beam_order(p, 20, 4, stop), "beam20");
        score(hollan::beam_order(p, 64, 4, stop), "beam64");
        score(hollan::beam_order(p, 128, 4, stop), "beam128");
    }
    for (const auto& [name, total] : totals)
        cout << left << setw(12) << name << " avg=" << double(total) / trials << '\n';
}
