#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "../sample_bots/hollan_bot.cpp"

using namespace std;

static int brute_best_cycles(const hollan::Problem& p) {
    vector<int> order(p.n);
    iota(order.begin(), order.end(), 0);
    int best = 0;
    do {
        vector<unsigned char> used(p.n, 0);
        bool valid = true;
        for (int k = 0; k < p.n; ++k) {
            int v = order[k];
            if (k) {
                bool touches = false;
                for (int u : p.nbr[v]) touches |= used[u];
                if (!touches) { valid = false; break; }
            }
            used[v] = 1;
        }
        if (valid) best = max(best, hollan::cycle_count(p, order));
    } while (next_permutation(order.begin(), order.end()));
    return best;
}

int main(int argc, char** argv) {
    int trials = argc > 1 ? stoi(argv[1]) : 30;
    mt19937_64 rng(0x455841435443484bULL);
    for (int t = 0; t < trials; ++t) {
        vector<int> values(9);
        iota(values.begin(), values.end(), 1);
        shuffle(values.begin(), values.end(), rng);
        vector<vector<int>> matrix(3, vector<int>(3));
        for (int i = 0; i < 9; ++i) matrix[i / 3][i % 3] = values[i];
        hollan::Problem p(matrix);
        int oracle = brute_best_cycles(p);
        hollan::Candidate exact = hollan::solve_for(p, 1.0);
        if (exact.cycles != oracle) {
            cerr << "mismatch trial=" << t << " exact=" << exact.cycles
                 << " oracle=" << oracle << '\n';
            return 2;
        }
    }
    cout << "exact_match=" << trials << '/' << trials << '\n';
}
