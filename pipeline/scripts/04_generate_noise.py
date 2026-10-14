#!/usr/bin/env python3
"""
Step 3.1 - Generate noisy variants of a query FASTA at several controlled levels.

Two noise types, matching the plan:
  substitution  - random base substitutions at rate r in {0.005, 0.01, 0.02, 0.05}
  truncation    - delete one random contiguous block of length f * genome_length,
                   f in {0.05, 0.10, 0.20}

Fixed seed per level for reproducibility. Each output FASTA keeps original
record ids so it can be joined back to metadata.csv unchanged.

Usage:
  python3 04_generate_noise.py --fasta data/T2.fasta --out-dir results/noise_T2 \
      --sub-rates 0.005,0.01,0.02,0.05 --trunc-fracs 0.05,0.10,0.20 --seed 42
"""
import argparse
import random
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

BASES = "ACGT"


def substitute(seq_str, rate, rng):
    seq = list(seq_str)
    n_sub = int(len(seq) * rate)
    positions = rng.sample(range(len(seq)), min(n_sub, len(seq)))
    for pos in positions:
        original = seq[pos].upper()
        choices = [b for b in BASES if b != original]
        seq[pos] = rng.choice(choices) if original in BASES else seq[pos]
    return "".join(seq)


def truncate(seq_str, frac, rng):
    n = len(seq_str)
    del_len = int(n * frac)
    if del_len <= 0 or del_len >= n:
        return seq_str
    start = rng.randint(0, n - del_len)
    return seq_str[:start] + seq_str[start + del_len:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sub-rates", default="0.005,0.01,0.02,0.05")
    ap.add_argument("--trunc-fracs", default="0.05,0.10,0.20")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(args.fasta, "fasta"))

    sub_rates = [float(x) for x in args.sub_rates.split(",") if x]
    trunc_fracs = [float(x) for x in args.trunc_fracs.split(",") if x]

    for rate in sub_rates:
        rng = random.Random(args.seed)
        out_path = out_dir / f"sub_{rate}.fasta"
        with open(out_path, "w") as f:
            for rec in records:
                noisy = substitute(str(rec.seq), rate, rng)
                rec2 = rec[:]
                rec2.seq = Seq(noisy)
                SeqIO.write(rec2, f, "fasta")
        print(f"[noise] wrote substitution level {rate} -> {out_path}")

    for frac in trunc_fracs:
        rng = random.Random(args.seed)
        out_path = out_dir / f"trunc_{frac}.fasta"
        with open(out_path, "w") as f:
            for rec in records:
                noisy = truncate(str(rec.seq), frac, rng)
                rec2 = rec[:]
                rec2.seq = Seq(noisy)
                SeqIO.write(rec2, f, "fasta")
        print(f"[noise] wrote truncation level {frac} -> {out_path}")

    # level 0 = clean copy, for a common x=0 baseline point on the robustness curve
    clean_path = out_dir / "level_0.fasta"
    with open(clean_path, "w") as f:
        SeqIO.write(records, f, "fasta")
    print(f"[noise] wrote clean baseline -> {clean_path}")


if __name__ == "__main__":
    main()
