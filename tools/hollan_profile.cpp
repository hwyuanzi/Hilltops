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

static uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static vector<vector<int>> matrix_for(int rows, int cols, uint64_t seed) {
    vector<int> a(rows * cols);
    iota(a.begin(), a.end(), 1);
    uint64_t state = seed;
    for (int i = rows * cols - 1; i > 0; --i)
        swap(a[i], a[splitmix64(state) % static_cast<uint64_t>(i + 1)]);
    vector<vector<int>> matrix(rows, vector<int>(cols));
    for (int i = 0; i < rows * cols; ++i) matrix[i / cols][i % cols] = a[i];
    return matrix;
}

int main(int argc, char** argv) {
    int trials = argc > 1 ? stoi(argv[1]) : 20;
    double seconds = argc > 2 ? stod(argv[2]) : 1.0;
    string path = argc > 3 ? argv[3] : "hollan_profile.csv";
    int rows = argc > 4 ? stoi(argv[4]) : 20;
    int cols = argc > 5 ? stoi(argv[5]) : rows;
    uint64_t base = 0xd1b54a32d192ed03ULL + 0x50000ULL;
    ofstream out(path);
    out << "index,seed,rows,cols,fallback,deterministic,post_local,post_beam,post_exact,final,"
           "det_wins,local_wins,beam_wins,exact_wins,repair_wins,restart_wins,destroy_wins,"
           "iterations,det_ms,local_ms,beam_ms,exact_ms,last_improvement_ms,total_ms\n";
    long long fallback = 0, deterministic = 0, local = 0, beam = 0, final = 0;
    long long det_wins = 0, local_wins = 0, beam_wins = 0, repair_wins = 0;
    long long restart_wins = 0, destroy_wins = 0, iterations = 0;
    double last_ms = 0.0, total_ms = 0.0;
    for (int i = 0; i < trials; ++i) {
        uint64_t seed = base + static_cast<uint64_t>(i);
        hollan::Problem p(matrix_for(rows, cols, seed));
        hollan::SearchStats s;
        auto start = chrono::steady_clock::now();
        hollan::Candidate candidate = hollan::solve_for(p, seconds, &s);
        double elapsed = 1000.0 * chrono::duration<double>(chrono::steady_clock::now() - start).count();
        if (!hollan::connected_order(p, candidate.order) || candidate.cycles != p.n - s.final_swaps) {
            cerr << "invalid profile trial " << i << '\n';
            return 2;
        }
        out << i << ',' << seed << ',' << rows << ',' << cols << ',' << s.fallback_swaps << ','
            << s.deterministic_swaps << ',' << s.post_local_swaps << ',' << s.post_beam_swaps << ','
            << s.post_exact_swaps << ',' << s.final_swaps << ',' << s.deterministic_wins << ','
            << s.local_wins << ',' << s.beam_wins << ',' << s.exact_wins << ','
            << s.random_repair_wins << ',' << s.random_restart_wins << ','
            << s.destroy_repair_wins << ',' << s.iterations << ',' << s.deterministic_ms << ','
            << s.local_ms << ',' << s.beam_ms << ',' << s.exact_ms << ','
            << s.last_improvement_ms << ',' << elapsed << '\n';
        fallback += s.fallback_swaps;
        deterministic += s.deterministic_swaps;
        local += s.post_local_swaps;
        beam += s.post_beam_swaps;
        final += s.final_swaps;
        det_wins += s.deterministic_wins;
        local_wins += s.local_wins;
        beam_wins += s.beam_wins;
        repair_wins += s.random_repair_wins;
        restart_wins += s.random_restart_wins;
        destroy_wins += s.destroy_repair_wins;
        iterations += s.iterations;
        last_ms += s.last_improvement_ms;
        total_ms += elapsed;
        cerr << "profile " << i + 1 << '/' << trials << " final=" << s.final_swaps << '\n';
    }
    auto avg = [trials](long long x) { return double(x) / trials; };
    cout << fixed << setprecision(3)
         << "trials=" << trials << " seconds=" << seconds
         << " fallback=" << avg(fallback) << " deterministic=" << avg(deterministic)
         << " local=" << avg(local) << " beam=" << avg(beam) << " final=" << avg(final) << '\n'
         << "wins/run det=" << avg(det_wins) << " local_cycles=" << avg(local_wins)
         << " beam=" << avg(beam_wins) << " random_repair=" << avg(repair_wins)
         << " restart=" << avg(restart_wins) << " destroy=" << avg(destroy_wins) << '\n'
         << "iterations/run=" << avg(iterations) << " last_improvement_ms=" << last_ms / trials
         << " total_ms=" << total_ms / trials << " output=" << path << '\n';
}
