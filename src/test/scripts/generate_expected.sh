#!/usr/bin/env bash
# Regenerates the expected alignment results in src/test/resources/expected
# by running the original C SSW library over the test dataset.
#
#   SSW_SRC=~/git/Complete-Striped-Smith-Waterman-Library/src src/test/scripts/generate_expected.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
SSW_SRC="${SSW_SRC:-$HOME/git/Complete-Striped-Smith-Waterman-Library/src}"
DATA="$ROOT/src/test/resources/dataset"
OUT="$ROOT/src/test/resources/expected"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

python3 -I "$ROOT/src/test/scripts/make_synthetic.py" "$DATA"

cc -O2 -Wall -I"$SSW_SRC" -o "$WORK/ssw_reference" "$ROOT/src/test/c/ssw_reference.c" "$SSW_SRC/ssw.c" -lm

cat "$DATA/uniprot_human.fasta" "$DATA/synthetic_queries.fasta" > "$WORK/queries.fasta"
cat "$DATA/uniprot_homologs.fasta" "$DATA/synthetic_targets.fasta" > "$WORK/targets.fasta"

mkdir -p "$OUT"
run() { # name type matrix gapO gapE flag targets queries
	local name=$1; shift
	"$WORK/ssw_reference" "$@" 2>/dev/null | gzip -9n > "$OUT/$name.tsv.gz"
	echo "$name: $(gzip -dc "$OUT/$name.tsv.gz" | wc -l) pairs"
}
run protein_blosum50_o3_e1 protein blosum50 3 1 1 "$WORK/targets.fasta" "$WORK/queries.fasta"
run protein_blosum50_o10_e2 protein blosum50 10 2 1 "$WORK/targets.fasta" "$WORK/queries.fasta"
run protein_blosum62_o11_e1 protein "$DATA/blosum62.txt" 11 1 1 "$WORK/targets.fasta" "$WORK/queries.fasta"
run dna_m2_x2_o3_e1 dna 2:2 3 1 1 "$DATA/dna_1k.fa" "$DATA/dna_54mer_hap1_1.100.fa"
