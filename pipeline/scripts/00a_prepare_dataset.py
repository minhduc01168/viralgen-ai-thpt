#!/usr/bin/env python3
"""
Step 0 (dataset prep) - Consolidate a "sars-cov2-time-forward-dataset/" style
directory tree into the data/T0.fasta ... T3.fasta + data/metadata.csv format
that scripts 00-06 expect.

Assumes the exact layout the user has:

  <root>/
    T0_train/
      sequences_<date>_<id>__<Lineage>_<PangoLineage>_<count>_<start>_<end>.fasta
      ...
    T1_test_pre_omicron/...
    T2_test_transition/...
    T3_test_post_omicron/...

Period is taken from the folder name (must start with T0/T1/T2/T3).
Lineage is parsed from the filename: the token right after "__".
Pango lineage / declared count / date range are also parsed and used only for
a sanity check (declared count in filename vs. actual sequences found).

Accession is taken as the first whitespace-separated token of the FASTA header,
with a trailing ".<version>" stripped (e.g. ">OQ123456.1 Severe acute..." ->
"OQ123456"). If your headers use a different scheme (e.g. GISAID EPI_ISL,
pipe-delimited), re-check with --dry-run before trusting the output.

Usage:
  python3 00_prepare_dataset.py --root data/sars-cov2-time-forward-dataset --outdir data
  python3 00_prepare_dataset.py --root data/sars-cov2-time-forward-dataset --outdir data --dry-run
"""
import argparse
import csv
import re
import sys
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

PERIOD_FOLDER_RE = re.compile(r"^(T[0-3])")
FNAME_RE = re.compile(
    r"^sequences_\d+_\d+__"
    r"(?P<lineage>[A-Za-z]+)_"
    r"(?P<pango>[A-Za-z0-9.]+)_"
    r"(?P<count>\d+)_"
    r"(?P<start>\d{1,2}[A-Za-z]{3}\d{4})_"
    r"(?P<end>\d{1,2}[A-Za-z]{3}\d{4})\.fasta?$"
)


def parse_filename(fname):
    m = FNAME_RE.match(fname)
    if not m:
        return None
    return m.groupdict()


def clean_accession(header_id):
    acc = header_id.split("|")[0]
    acc = acc.split(".")[0]  # drop version suffix
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="path to sars-cov2-time-forward-dataset/")
    ap.add_argument("--outdir", default="data")
    ap.add_argument("--dry-run", action="store_true", help="parse & report only, write nothing")
    args = ap.parse_args()

    root = Path(args.root)
    if not root.is_dir():
        sys.exit(f"[ERROR] root not found: {root}")

    outdir = Path(args.outdir)
    if not args.dry_run:
        outdir.mkdir(parents=True, exist_ok=True)

    period_dirs = sorted(p for p in root.iterdir() if p.is_dir())
    if not period_dirs:
        sys.exit(f"[ERROR] no subdirectories found under {root}")

    metadata_rows = []
    seen_accessions = set()
    duplicate_accessions = []
    unparsed_files = []
    period_seq_counts = defaultdict(int)
    period_writers = {}
    period_out_paths = {}

    for pdir in period_dirs:
        m = PERIOD_FOLDER_RE.match(pdir.name)
        if not m:
            print(f"[skip] folder does not start with T0/T1/T2/T3: {pdir.name}")
            continue
        period = m.group(1)

        fasta_files = sorted(pdir.glob("*.fasta")) + sorted(pdir.glob("*.fa"))
        if not fasta_files:
            print(f"[warn] no .fasta files found in {pdir}")
            continue

        if not args.dry_run and period not in period_writers:
            out_path = outdir / f"{period}.fasta"
            period_writers[period] = open(out_path, "w")
            period_out_paths[period] = out_path

        for fpath in fasta_files:
            info = parse_filename(fpath.name)
            if info is None:
                unparsed_files.append(str(fpath))
                print(f"[WARN] could not parse lineage from filename, skipping: {fpath.name}")
                continue

            lineage = info["lineage"]
            declared_count = int(info["count"])

            records = list(SeqIO.parse(fpath, "fasta"))
            actual_count = len(records)
            status = "OK" if actual_count == declared_count else "MISMATCH"
            print(f"[{period}/{lineage}] {fpath.name}: declared={declared_count} actual={actual_count} [{status}]")

            for rec in records:
                acc = clean_accession(rec.id)
                if acc in seen_accessions:
                    duplicate_accessions.append(acc)
                    continue
                seen_accessions.add(acc)
                metadata_rows.append({
                    "accession": acc,
                    "lineage": lineage,
                    "period": period,
                    "collection_date": f"{info['start']}-{info['end']}",
                })
                period_seq_counts[period] += 1
                if not args.dry_run:
                    SeqIO.write(rec, period_writers[period], "fasta")

    for w in period_writers.values():
        w.close()

    print("\n=== Summary ===")
    for period in sorted(period_seq_counts):
        print(f"  {period}: {period_seq_counts[period]} sequences written"
              + ("" if args.dry_run else f" -> {period_out_paths[period]}"))
    if duplicate_accessions:
        print(f"[WARNING] {len(duplicate_accessions)} duplicate accessions were skipped "
              f"(kept first occurrence only): {duplicate_accessions[:10]}"
              + (" ..." if len(duplicate_accessions) > 10 else ""))
    if unparsed_files:
        print(f"[WARNING] {len(unparsed_files)} file(s) had a name that didn't match the "
              f"expected pattern and were SKIPPED — check these manually:")
        for f in unparsed_files:
            print(f"   - {f}")

    if args.dry_run:
        print("\n[dry-run] no files written. Re-run without --dry-run once the summary above looks correct.")
        return

    meta_path = outdir / "metadata.csv"
    with open(meta_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "lineage", "period", "collection_date"])
        w.writeheader()
        w.writerows(metadata_rows)
    print(f"\n[done] wrote {len(metadata_rows)} rows to {meta_path}")
    print("[done] period FASTA files written to:", ", ".join(str(p) for p in period_out_paths.values()))


if __name__ == "__main__":
    main()
