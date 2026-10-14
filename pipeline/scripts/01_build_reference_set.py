#!/usr/bin/env python3
"""
Step 2.2 - Build a small, representative reference set from T0.

Input expected:
  --fasta   data/T0.fasta            (all T0 genomes, one record per accession)
  --meta    data/metadata.csv        (columns: accession,lineage,period,collection_date)
            period must contain 'T0' for training-period rows.
Output:
  results/reference.fasta            (n_per_class representative sequences per lineage)
  results/reference_metadata.csv

Selection strategy: stratified random sampling per lineage, fixed seed (matches
the paper's principle of fixing everything on T0 before touching later periods).
Using random sampling (not centroid-based) keeps this auditable/reproducible with
a single seed; switch to k-mer-centroid selection later if reviewers want it.
"""
import argparse
import random
import csv
from pathlib import Path
from Bio import SeqIO


def load_metadata(meta_path):
    meta = {}
    with open(meta_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            meta[row["accession"]] = row
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="T0 FASTA file")
    ap.add_argument("--meta", required=True, help="metadata CSV (accession,lineage,period,collection_date)")
    ap.add_argument("--n-per-class", type=int, default=30, help="reference sequences per lineage")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    random.seed(args.seed)
    meta = load_metadata(args.meta)

    # group T0 accessions by lineage
    by_lineage = {}
    for acc, row in meta.items():
        if row["period"] != "T0":
            continue
        by_lineage.setdefault(row["lineage"], []).append(acc)

    chosen = set()
    for lineage, accs in by_lineage.items():
        k = min(args.n_per_class, len(accs))
        picked = random.sample(accs, k)
        chosen.update(picked)
        print(f"[reference] {lineage}: picked {k}/{len(accs)}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ref_fasta = outdir / "reference.fasta"
    ref_meta = outdir / "reference_metadata.csv"

    n_written = 0
    with open(ref_fasta, "w") as fout, open(ref_meta, "w", newline="") as mout:
        writer = csv.writer(mout)
        writer.writerow(["accession", "lineage"])
        for rec in SeqIO.parse(args.fasta, "fasta"):
            acc = rec.id.split(".")[0].split("|")[0]
            if acc in chosen:
                SeqIO.write(rec, fout, "fasta")
                writer.writerow([acc, meta[acc]["lineage"]])
                n_written += 1

    print(f"[reference] wrote {n_written} sequences to {ref_fasta}")
    if n_written < len(chosen):
        print(f"[WARNING] {len(chosen) - n_written} chosen accessions were not found "
              f"in the FASTA file — check accession formatting (version suffixes, etc.)")


if __name__ == "__main__":
    main()
