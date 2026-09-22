#!/usr/bin/env python3
import csv
import statistics
import sys


def load(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return rows


def percentile(values, p):
    values = sorted(values)
    if not values:
        return float("nan")
    return values[min(len(values) - 1, int(p * len(values)))]


def summarize(label, rows):
    swaps = [int(r["swaps"]) for r in rows]
    times = [float(r["runtime_ms"]) for r in rows]
    print(
        f"{label}: n={len(rows)} mean={statistics.mean(swaps):.4f} "
        f"median={statistics.median(swaps):.2f} p75={percentile(swaps,.75)} "
        f"p90={percentile(swaps,.90)} p95={percentile(swaps,.95)} "
        f"min={min(swaps)} max={max(swaps)} runtime_mean_ms={statistics.mean(times):.3f} "
        f"runtime_max_ms={max(times):.3f} valid={sum(int(r['valid']) for r in rows)}/{len(rows)}"
    )


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: hollan_compare.py BASE.csv CANDIDATE.csv")
    base, cand = load(sys.argv[1]), load(sys.argv[2])
    if len(base) != len(cand):
        raise SystemExit("row count differs")
    for a, b in zip(base, cand):
        key_a = (a["index"], a["seed"], a["rows"], a["cols"], a["family"])
        key_b = (b["index"], b["seed"], b["rows"], b["cols"], b["family"])
        if key_a != key_b:
            raise SystemExit(f"corpus mismatch: {key_a} != {key_b}")
    summarize("baseline", base)
    summarize("candidate", cand)
    delta = [int(b["swaps"]) - int(a["swaps"]) for a, b in zip(base, cand)]
    wins = sum(d < 0 for d in delta)
    ties = sum(d == 0 for d in delta)
    losses = sum(d > 0 for d in delta)
    print(
        f"paired: mean_delta={statistics.mean(delta):+.4f} wins/ties/losses={wins}/{ties}/{losses} "
        f"win2={sum(d <= -2 for d in delta)} loss2={sum(d >= 2 for d in delta)} "
        f"best_delta={min(delta):+d} worst_regression={max(delta):+d}"
    )
    for size in sorted({(int(r["rows"]), int(r["cols"])) for r in base}):
        idx = [i for i, r in enumerate(base) if (int(r["rows"]), int(r["cols"])) == size]
        ds = [delta[i] for i in idx]
        print(
            f"  {size[0]}x{size[1]} n={len(idx)} delta={statistics.mean(ds):+.4f} "
            f"W/T/L={sum(d<0 for d in ds)}/{sum(d==0 for d in ds)}/{sum(d>0 for d in ds)} "
            f"loss2={sum(d>=2 for d in ds)}"
        )


if __name__ == "__main__":
    main()
