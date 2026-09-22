# Hollan bot experiment log

Starting solver: `7a157e8c54874bb15e353c833f1f5d6e4ebd7001`.

All paired comparisons use `hollan_corpus_bench.cpp` and identical deterministic
matrices. The full frozen development corpus contains the requested 2,400 square
random cases plus 600 thin/rectangular cases and 25 structured patterns. The
held-out corpus uses a disjoint seed range and is not used for parameter choices.

| Experiment | Hypothesis | Change | Corpus | Result | Decision |
|---|---|---|---|---|---|
| Frozen baseline | Establish reproducible reference | None | Full: 3,025 cases, 20 ms | 3,025/3,025 valid; mean 22.3025, median 14, p90 60, p95 107; max runtime 32.523 ms | Reference |
| Skip large beam (timed pilot) | 20x20 beam's 12% time slice yielded no wins in attribution | Run beam only for `N <= 225` | Full, 20 ms | On 20x20: 2 wins/203 ties/0 losses; unrelated groups showed timer noise | Retest with fixed work |
| Replace full restarts | General randomized restarts won 0/30 profiled 20x20 runs; prefix rebuild won on 23/30 | Make all randomized greedy constructions prefix destroy/rebuild for `N >= 8` | 20x20 train: 50 cases, 512 fixed iterations | Mean +0.12 swaps; 15/20/15 W/T/L; best -9 but worst +10 | Reject; poor tail robustness |
| Complementary repair search | A repair rule should optimize incremental cycles beyond the whole deterministic portfolio | Random-search 1,000 parameter sets on 100 train matrices; retain top 3 rules without changing random bases; apply for `N >= 60` | 50 train + 50 held-out 20x20, 512 fixed iterations; full 3,025-case timed corpus | Held-out 20x20: -12.56 mean swaps, 49/1/0 W/T/L, no regression. Full timed: -1.954 mean overall, -14.663 on 20x20 | Keep |
