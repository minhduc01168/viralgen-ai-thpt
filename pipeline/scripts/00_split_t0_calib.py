#!/usr/bin/env python3
"""
Helper - build T0_calib.fasta = T0 minus the reference-set accessions.

Why: the anomaly threshold must be fit on IND data that is NOT the reference
set itself (otherwise every reference sequence trivially matches itself at
100% identity / distance 0, giving a degenerate threshold). This mirrors how
the paper fits theta on T0 Mahalanobis distances without leaking test data.

Usage:
  python3 00_split_t0_calib.py --fasta data/T0.fasta \
      --ref-meta results/reference_metadata.csv \
      --out results/T0_calib.fasta
"""
import argparse
import csv
from pathlib import Path
from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--ref-meta", required=True)
    ap.add_argument("--out", default="results/T0_calib.fasta")
    args = ap.parse_args()

    with open(args.ref_meta) as f:
        ref_accs = {row["accession"] for row in csv.DictReader(f)}

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(args.out, "w") as fout:
        for rec in SeqIO.parse(args.fasta, "fasta"):
            acc = rec.id.split(".")[0].split("|")[0]
            if acc not in ref_accs:
                SeqIO.write(rec, fout, "fasta")
                n += 1
    print(f"[calib] wrote {n} T0 sequences (excluding {len(ref_accs)} reference accessions) to {args.out}")


if __name__ == "__main__":
    main()
