#!/usr/bin/env python3
"""
Step 3.3 - Aggregate robustness-to-noise results into a CSV and Figure 7.

For each noise level, computes the rejection accuracy of each baseline against
the ground-truth label (known lineage should NOT be rejected; Omicron SHOULD
be rejected), using the threshold already fit in Step 2 (table7_comparison.csv
holds the fitted threshold per method).

IMPORTANT: this script produces the BLAST/Mash curves only. To finish Figure 7
as planned (two-line comparison), re-run your own encoder + Mahalanobis
pipeline on the exact same noisy FASTA files in results/noise_<PERIOD>/ and
append its accuracy numbers to results/robustness_summary.csv under
method="encoder" before re-running the plotting step below.

Usage:
  python3 06_aggregate_robustness.py --meta data/metadata.csv --period T2 --method blast
  python3 06_aggregate_robustness.py --meta data/metadata.csv --period T2 --method mash
  python3 06_aggregate_robustness.py --plot-only   # after all methods (incl. encoder) are in the CSV
"""
import argparse
import csv
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_metadata(meta_path):
    meta = {}
    with open(meta_path) as f:
        for row in csv.DictReader(f):
            meta[row["accession"]] = row
    return meta


def best_hit(tsv_path, method):
    best = {}
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if method == "blast":
                if len(parts) < 3:
                    continue
                q, pident = parts[0].split(".")[0].split("|")[0], float(parts[2])
                if q not in best or pident > best[q]:
                    best[q] = pident
            else:
                if len(parts) < 3:
                    continue
                q, dist = parts[1].split(".")[0].split("|")[0], float(parts[2])
                if q not in best or dist < best[q]:
                    best[q] = dist
    return best


def novelty(best_hits, method):
    if method == "blast":
        return {k: 100.0 - v for k, v in best_hits.items()}
    return dict(best_hits)


def level_sort_key(level):
    m = re.match(r"(sub|trunc|level)_?([\d.]*)", level)
    kind, val = m.group(1), m.group(2)
    order = {"level": -1, "sub": 0, "trunc": 1}
    return (order.get(kind, 2), float(val) if val else 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta")
    ap.add_argument("--period")
    ap.add_argument("--method", choices=["blast", "mash"])
    ap.add_argument("--resdir", default="results")
    ap.add_argument("--plot-only", action="store_true")
    args = ap.parse_args()

    summary_path = Path(args.resdir) / "robustness_summary.csv"

    if not args.plot_only:
        meta = load_metadata(args.meta)
        threshold = None
        with open(Path(args.resdir) / "table7_comparison.csv") as f:
            for row in csv.DictReader(f):
                if row["method"] == args.method:
                    threshold = float(row["threshold"])
        if threshold is None:
            raise SystemExit(f"No fitted threshold for method={args.method} in table7_comparison.csv "
                              f"— run 03_compute_metrics.py first.")

        robdir = Path(args.resdir) / f"robustness_{args.period}"
        levels = sorted(
            {p.stem.split(f"{args.method}_")[1] for p in robdir.glob(f"{args.method}_*.tsv")},
            key=level_sort_key,
        )

        rows = []
        for level in levels:
            hits = best_hit(robdir / f"{args.method}_{level}.tsv", args.method)
            scores = novelty(hits, args.method)

            correct = 0
            total = 0
            for acc, row in meta.items():
                if row["period"] != args.period:
                    continue
                total += 1
                score = scores.get(acc, np.inf)
                predicted_unknown = score > threshold
                is_omicron = row["lineage"].lower() == "omicron"
                if predicted_unknown == is_omicron:
                    correct += 1
            acc_rate = correct / total if total else float("nan")
            rows.append((level, acc_rate, total))
            print(f"[{args.method}] {level}: accuracy={acc_rate:.4f} (n={total})")

        write_header = not summary_path.exists()
        with open(summary_path, "a", newline="") as f:
            w = csv.writer(f)
            if write_header:
                w.writerow(["method", "period", "level", "accuracy", "n"])
            for level, acc_rate, total in rows:
                w.writerow([args.method, args.period, level, round(acc_rate, 4), total])
        print(f"[{args.method}] appended to {summary_path}")

    # --- plotting: one line per method found in the summary CSV ---
    if not summary_path.exists():
        print("[plot] no summary CSV yet, skipping plot.")
        return

    by_method = {}
    with open(summary_path) as f:
        for row in csv.DictReader(f):
            by_method.setdefault(row["method"], []).append(
                (level_sort_key(row["level"]), row["level"], float(row["accuracy"]))
            )

    plt.figure(figsize=(6, 4))
    for method, points in by_method.items():
        points.sort(key=lambda x: x[0])
        labels = [p[1] for p in points]
        accs = [p[2] for p in points]
        plt.plot(labels, accs, marker="o", label=method)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Rejection accuracy")
    plt.xlabel("Noise level")
    plt.title("Robustness to noise: encoder vs. alignment-based baselines")
    plt.legend()
    plt.tight_layout()
    fig_path = Path(args.resdir).parent / "figures" / "figure7_robustness_T2.png"
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig_path, dpi=200)
    print(f"[plot] wrote {fig_path}")


if __name__ == "__main__":
    main()
