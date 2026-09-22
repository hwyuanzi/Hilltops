# Hollan bot experiment log

Starting solver: `7a157e8c54874bb15e353c833f1f5d6e4ebd7001`.

All paired comparisons use `hollan_corpus_bench.cpp` and identical deterministic
matrices. The full frozen development corpus contains the requested 2,400 square
random cases plus 600 thin/rectangular cases and 25 structured patterns. The
held-out corpus uses a disjoint seed range and is not used for parameter choices.

| Experiment | Hypothesis | Change | Corpus | Result | Decision |
|---|---|---|---|---|---|
| Frozen baseline | Establish reproducible reference | None | Full: 3,025 cases, 20 ms | 3,025/3,025 valid; mean 22.3025, median 14, p90 60, p95 107; max runtime 32.523 ms | Reference |
