# Hilltops competition budget experiments

Baseline: commit `7baa707`. The competition runner is copied from Emil's
`21e827d` (`TIMEOUT_SECONDS = 120`, a new `Bot()` per entrant, Julia `-t auto`).
All bot experiments used the frozen matrices in `hollan_corpus_bench.cpp`,
compiled with `c++ -std=c++17 -O3`, and checked the returned swaps with the
independent Hilltop evaluator. Lower swap counts are better. Paired timings are
wall-clock solver budgets on this machine, not performance guarantees for the
architect's machine. The submitted `get_swaps` budget is 112 seconds, leaving
eight seconds under the runner's hard limit; the 120-second row below tests the
search method directly, outside the submission wrapper.

## Paired 20x20 validation budgets

| Seconds | Boards | Baseline mean swaps | Final mean swaps | Final W/T/L |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 10 | 84.30 | 83.70 | 3/5/2 |
| 10 | 6 | 81.50 | 80.67 | 3/1/2 |
| 30 | 3 | 81.33 | 80.00 | 2/1/0 |
| 60 | 2 | 82.00 | 79.50 | 2/0/0 |
| 120 | 2 | 81.50 | 79.50 | 2/0/0 |

The rows use the first *n* boards of the same validation corpus. These are
small samples at longer budgets; means across rows should not be interpreted
as a single board's improvement curve.

## Held-out 20x20 boards

| Seconds | Boards | Baseline mean swaps | Final mean swaps | Final W/T/L |
| ---: | ---: | ---: | ---: | ---: |
| 10 | 8 | 82.00 | 81.25 | 6/0/2 |
| 30 | 3 | 82.00 | 82.00 | 1/1/1 |
| 60 | 2 | 83.00 | 82.50 | 1/1/0 |
| 120 | 2 | 82.50 | 81.50 | 1/1/0 |

## One-change experiments

- **Keep: randomized repair-beam choice jitter.** Add uniform score jitter of
  magnitude 5 only to the anytime randomized repair beams; leave initial
  deterministic beams unchanged. It improved held-out 10-second quality and
  the tested 30/60/120-second validation means.
- **Reject: 50% late beam allocation.** Replace greedy restart slots with
  repair beams after 10 seconds. On three 30-second validation boards, the
  mean was unchanged, with one win, one tie, and one loss.
- **Reject: cached cycle membership.** Mark the cycle containing the next
  position once per beam state instead of walking it for each child. Ten
  fixed-256-iteration boards tied exactly, while mean runtime rose from
  381.69 to 383.60 ms.
- **Reject: annealing neighborhood.** Spend 2 seconds on the baseline and
  8 seconds on connected-order transposition annealing. On six 10-second
  validation boards it averaged 82.67 swaps versus 81.50 for uninterrupted
  baseline search (zero wins, two ties, four losses).
- **Reject: success-conditioned beam weights.** Mutate weights around the
  last winning beam on half of later beams. It improved eight held-out
  10-second boards against the baseline (82.00 to 81.125), but the retained
  jitter variant averaged 81.25 there and was better on three 30-second
  validation boards (80.00 versus 81.00). The extra mechanism had no clear
  incremental benefit.

Final correctness: 3,025/3,025 full frozen-corpus cases valid at 20 ms using
independent swap replay, connected-target, cycle-count, and worst-distance
checks. The short-budget score is not the selection target for this competition
run. One end-to-end 20x20 invocation through Emil's runner compiled and ran
the C++ bot, returned 80 swaps, passed an independent Hilltop check, and took
113.688 seconds total, below the runner's 120-second timeout.

## Post-checkpoint search, commit `7065b6f`

The checkpoint was pushed to `origin/main` before these trials. Each timed
variant was compared against that frozen executable on identical held-out
20x20 matrices; `hollan_compare.py` checked corpus identity and independent
swap validity. The held-out set has now been used to select variants, so it is
no longer an untouched final test set. Negative deltas are fewer swaps.

| Single conceptual change | 1 s, 50 boards: delta and W/T/L | 10 s, 8 boards: delta and W/T/L | Decision |
| --- | ---: | ---: | --- |
| Seed repair beams from the incumbent and branch at connected cycle splits | -0.02; 22/12/16 | -0.375; 4/2/2 | Reject: mixed and too small |
| Locally improve promising random-beam candidates before incumbent comparison | -0.10; 20/11/19 | -0.125; 3/3/2 | Reject: mixed and too small |
| Raise random-beam choice noise from 5 to 25 after 128 stagnant iterations | +1.36; 0/22/28 | Stopped after clear loss | Reject |
| Broaden repair-rule mix after 128 stagnant iterations | -0.08; 18/17/15 | -0.25; 2/3/3 | Reject: mixed and too small |
| Bias cycle-aware transposition attempts toward positions within two board widths | +0.32; 17/14/19 | 0.00; 4/0/4 | Reject |
| Locally improve all near-top completed repair-beam paths | +0.06; 15/21/14 | +0.125; 1/5/2 | Reject |

For a true compound neighborhood, an exhaustive probe paired transpositions
from distinct permutation cycles on 12 held-out one-second incumbents. It
searched for two-swap cycle gains where the first swap was disconnected but
the final target was connected. No such gain appeared on any of the 12 boards;
the operator was not added to the timed solver.

All search-method variants above were reverted. The only submitted-code change
after the checkpoint reduces the search budget from 112 to 108 seconds, giving
the 120-second runner twelve seconds for compilation, I/O, scheduling, and
submission. Its score tradeoff has not been established by a large paired
long-budget sample; the change is for timing safety.

Final checks for the 108-second file: 3,025/3,025 frozen-corpus boards valid
at a 20 ms solver budget. A separate 20x20 call through `Bot.get_swaps` compiled
the single C++ file, returned 80 swaps, passed an independent Hilltop replay,
and took 109.826 seconds end to end on this machine. The paired benchmark CSVs
and final correctness CSV are saved under `results/hollan_evolution/`.
