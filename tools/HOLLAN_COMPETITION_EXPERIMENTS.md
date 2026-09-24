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
