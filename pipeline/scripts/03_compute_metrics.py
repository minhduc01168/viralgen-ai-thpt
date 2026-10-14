#!/usr/bin/env python3
"""
Step 2.3-2.5 - Fit threshold on T0_calib, classify T1/T2/T3, compute:
  - per-lineage-per-period rejection rate  (mirrors paper Table 6)
  - AUROC / FPR95 with T1 = IND, pooled Omicron(T2+T3) = OOD (mirrors Section 3.6 / Table 4)
  - runtime summary (mirrors the "Runtime" column of the new Table 7)

Run AFTER 02_run_baselines.sh has produced blast_/mash_ TSVs for T0_calib, T1, T2, T3.

Usage:
  python3 03_compute_metrics.py --meta data/metadata.csv --method blast
  python3 03_compute_metrics.py --meta data/metadata.csv --method mash
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score


def load_metadata(meta_path):
    meta = {}
    with open(meta_path) as f:
        for row in csv.DictReader(f):
            meta[row["accession"]] = row
    return meta


def best_hit_blast(tsv_path):
    """qseqid -> best (highest) pident among hits."""
    best = {}
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            q, s, pident = parts[0], parts[1], float(parts[2])
            q = q.split(".")[0].split("|")[0]
            if q not in best or pident > best[q]:
                best[q] = pident
    return best  # higher = more similar (score, not distance)


def best_hit_mash(tsv_path):
    """query -> best (lowest) mash distance among hits."""
    best = {}
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            ref, q, dist = parts[0], parts[1], float(parts[2])
            q = q.split(".")[0].split("|")[0]
            if q not in best or dist < best[q]:
                best[q] = dist
    return best  # lower = more similar (a true distance)


def to_novelty_score(best_hits, method):
    """Return dict accession -> novelty score, HIGHER = more anomalous, for both methods."""
    if method == "blast":
        # identity is a similarity; novelty = 100 - pident (missing hit => 100, maximally novel)
        return {k: (100.0 - v) for k, v in best_hits.items()}
    else:
        return dict(best_hits)  # mash distance already higher = more novel


def get_runtime(timing_path):
    secs = {}
    if not Path(timing_path).exists():
        return secs
    with open(timing_path) as f:
        for line in f:
            k, v = line.strip().split("=")
            secs[k] = float(v)
    return secs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--method", choices=["blast", "mash"], required=True)
    ap.add_argument("--resdir", default="results")
    ap.add_argument("--percentile", type=float, default=95.0,
                     help="threshold = this percentile of T0_calib novelty scores (matches theta at 95th pct in the paper)")
    args = ap.parse_args()

    meta = load_metadata(args.meta)
    resdir = Path(args.resdir)

    def hits_for(period):
        path = resdir / f"{args.method}_{period}.tsv"
        parser = best_hit_blast if args.method == "blast" else best_hit_mash
        return parser(path)

    # --- fit threshold on T0_calib ---
    calib_hits = hits_for("T0_calib")
    calib_scores = to_novelty_score(calib_hits, args.method)
    threshold = np.percentile(list(calib_scores.values()), args.percentile)
    print(f"[{args.method}] threshold (novelty score, {args.percentile}th pct of T0_calib) = {threshold:.4f}")
    print(f"[{args.method}] T0_calib n={len(calib_scores)}, missing (no hit)={0}")

    # --- classify T1/T2/T3, mirror Table 6 ---
    summary_rows = []
    period_scores = {}
    for period in ["T1", "T2", "T3"]:
        hits = hits_for(period)
        scores = to_novelty_score(hits, args.method)
        period_scores[period] = scores

        by_lineage = {}
        for acc, row in meta.items():
            if row["period"] != period:
                continue
            lineage = row["lineage"]
            score = scores.get(acc, np.inf)  # no BLAST/Mash hit at all => maximally novel
            by_lineage.setdefault(lineage, []).append(score)

        for lineage, scores_list in sorted(by_lineage.items()):
            arr = np.array(scores_list)
            rejected = (arr > threshold).sum()
            n = len(arr)
            rate = rejected / n if n else float("nan")
            summary_rows.append({
                "period": period, "lineage": lineage, "n": n,
                "rejected": int(rejected), "rejection_rate": round(rate, 4),
                "mean_novelty_score": round(float(arr.mean()), 4),
            })
            print(f"  {period} {lineage}: n={n} rejection_rate={rate:.3f}")

    out_csv = resdir / f"table6_{args.method}.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["period", "lineage", "n", "rejected", "rejection_rate", "mean_novelty_score"])
        w.writeheader()
        w.writerows(summary_rows)
    print(f"[{args.method}] wrote {out_csv}")

    # --- AUROC / FPR95: IND = T1, OOD = pooled Omicron from T2+T3 (mirrors Section 3.6) ---
    y_true, y_score = [], []
    for acc, row in meta.items():
        if row["period"] == "T1":
            y_true.append(0)
            y_score.append(period_scores["T1"].get(acc, np.inf))
        elif row["period"] in ("T2", "T3") and row["lineage"].lower() == "omicron":
            y_true.append(1)
            y_score.append(period_scores[row["period"]].get(acc, np.inf))

    y_true = np.array(y_true)
    y_score = np.array(y_score, dtype=float)
    finite_max = y_score[np.isfinite(y_score)].max() if np.isfinite(y_score).any() else 1.0
    y_score = np.where(np.isfinite(y_score), y_score, finite_max * 10)  # push "no hit" beyond max

    auroc = roc_auc_score(y_true, y_score)
    # FPR95: FPR when TPR (recall on OOD) = 95%
    order = np.argsort(-y_score)
    y_true_sorted = y_true[order]
    tp_cum = np.cumsum(y_true_sorted)
    fp_cum = np.cumsum(1 - y_true_sorted)
    n_pos = y_true.sum()
    n_neg = (1 - y_true).sum()
    tpr = tp_cum / n_pos
    fpr = fp_cum / n_neg
    idx = np.searchsorted(tpr, 0.95)
    idx = min(idx, len(fpr) - 1)
    fpr95 = fpr[idx]

    print(f"[{args.method}] AUROC (T1 IND vs pooled Omicron OOD) = {auroc:.4f}")
    print(f"[{args.method}] FPR95 = {fpr95:.4f}")

    # --- runtime summary across periods ---
    total_secs = 0.0
    for period in ["T1", "T2", "T3"]:
        t = get_runtime(resdir / f"timing_{period}.txt")
        total_secs += sum(t.values())
    print(f"[{args.method}] total wall-clock over T1+T2+T3 = {total_secs:.2f}s")

    # --- append to master comparison table (Table 7) ---
    table7_path = resdir / "table7_comparison.csv"
    write_header = not table7_path.exists()
    with open(table7_path, "a", newline="") as f:
        w = csv.writer(f)
        if write_header:
            w.writerow(["method", "auroc", "fpr95", "runtime_seconds_T1_T3", "threshold"])
        w.writerow([args.method, round(auroc, 4), round(fpr95, 4), round(total_secs, 2), round(threshold, 4)])
    print(f"[{args.method}] appended row to {table7_path}")


if __name__ == "__main__":
    main()
