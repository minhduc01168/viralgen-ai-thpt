#!/usr/bin/env bash
# Step 2.3 - Run BLAST and Mash baselines: reference (T0) vs each query period.
#
# Usage:
#   bash 02_run_baselines.sh T0 results/reference.fasta   # self-distance, for threshold fitting
#   bash 02_run_baselines.sh T1 data/T1.fasta
#   bash 02_run_baselines.sh T2 data/T2.fasta
#   bash 02_run_baselines.sh T3 data/T3.fasta
#
# Outputs (all under results/):
#   blast_<PERIOD>.tsv   - outfmt 6: qseqid sseqid pident length evalue bitscore
#   mash_<PERIOD>.tsv    - mash dist output: ref query dist p-value shared-hashes
#   timing_<PERIOD>.txt  - wall-clock seconds for each tool
set -euo pipefail

PERIOD=$1
QUERY_FASTA=$2
REF_FASTA=${3:-results/reference.fasta}
OUTDIR=results
mkdir -p "$OUTDIR"

echo "[baseline] period=$PERIOD query=$QUERY_FASTA ref=$REF_FASTA"

# --- BLAST ---
if [ ! -f "$OUTDIR/refdb.nin" ] && [ ! -f "$OUTDIR/refdb.ndb" ]; then
    makeblastdb -in "$REF_FASTA" -dbtype nucl -out "$OUTDIR/refdb"
fi

START=$(date +%s.%N)
blastn -query "$QUERY_FASTA" -db "$OUTDIR/refdb" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -max_target_seqs 1 -num_threads "$(nproc)" \
    -out "$OUTDIR/blast_${PERIOD}.tsv"
END=$(date +%s.%N)
echo "blast_seconds=$(echo "$END - $START" | bc)" > "$OUTDIR/timing_${PERIOD}.txt"

# --- Mash ---
# -i = sketch each sequence in the multi-FASTA individually (default sketches
# the whole file as one combined sketch, which is wrong for per-genome queries)
if [ ! -f "$OUTDIR/refdb.msh" ]; then
    mash sketch -i -o "$OUTDIR/refdb" "$REF_FASTA"
fi
mash sketch -i -o "$OUTDIR/query_${PERIOD}" "$QUERY_FASTA"

START=$(date +%s.%N)
mash dist "$OUTDIR/refdb.msh" "$OUTDIR/query_${PERIOD}.msh" > "$OUTDIR/mash_${PERIOD}.tsv"
END=$(date +%s.%N)
echo "mash_seconds=$(echo "$END - $START" | bc)" >> "$OUTDIR/timing_${PERIOD}.txt"

echo "[baseline] done -> $OUTDIR/blast_${PERIOD}.tsv , $OUTDIR/mash_${PERIOD}.tsv"
cat "$OUTDIR/timing_${PERIOD}.txt"
