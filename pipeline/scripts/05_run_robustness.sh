#!/usr/bin/env bash
# Step 3.2 - Run BLAST/Mash baselines across all noise levels for a given period.
#
# Usage: bash 05_run_robustness.sh T2 data/T2.fasta
#
# Produces, per level (e.g. sub_0.01, trunc_0.10, level_0):
#   results/robustness_<PERIOD>/blast_<LEVEL>.tsv
#   results/robustness_<PERIOD>/mash_<LEVEL>.tsv
set -euo pipefail

PERIOD=$1
QUERY_FASTA=$2
REF_FASTA=results/reference.fasta
NOISE_DIR="results/noise_${PERIOD}"
OUTDIR="results/robustness_${PERIOD}"
mkdir -p "$OUTDIR"

# 1. generate noisy variants
python3 scripts/04_generate_noise.py --fasta "$QUERY_FASTA" --out-dir "$NOISE_DIR"

# 2. run both baselines at every noise level
for LEVEL_FASTA in "$NOISE_DIR"/*.fasta; do
    LEVEL=$(basename "$LEVEL_FASTA" .fasta)
    echo "=== level: $LEVEL ==="

    blastn -query "$LEVEL_FASTA" -db results/refdb \
        -outfmt "6 qseqid sseqid pident length evalue bitscore" \
        -max_target_seqs 1 -num_threads "$(nproc)" \
        -out "$OUTDIR/blast_${LEVEL}.tsv"

    mash sketch -i -o "$OUTDIR/query_${LEVEL}" "$LEVEL_FASTA" 2>/dev/null
    mash dist results/refdb.msh "$OUTDIR/query_${LEVEL}.msh" > "$OUTDIR/mash_${LEVEL}.tsv"
done

echo "[robustness] done -> $OUTDIR"
