#!/usr/bin/env bash
# Checks the native library (target/native/libcssw4j.so, built with `mvn -Pnative package`) against the reference
# results of the original C SSW library in src/test/resources/expected.
#
#   src/test/scripts/test_native.sh
#
# Columns compared: query, target, score1, begin/end positions, flag and cigar. The second best score and its
# position (columns 4 and 9) depend on the vector width and are not compared, as in CReferenceTest.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
NATIVE="$ROOT/target/native"
DATA="$ROOT/src/test/resources/dataset"
EXPECTED="$ROOT/src/test/resources/expected"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [[ ! -f "$NATIVE/libcssw4j.so" ]]; then
	echo "error: $NATIVE/libcssw4j.so not found; build it with: GRAALVM_HOME=... mvn -Pnative package" >&2
	exit 1
fi

cc -O2 -Wall -Wno-unused-function -Wno-unused-variable -I"$NATIVE" \
	-o "$WORK/ssw_native_check" "$ROOT/src/test/native/ssw_native_check.c" \
	-L"$NATIVE" -l:libcssw4j.so -Wl,-rpath,"$NATIVE" -lm

COMPARE='1,2,3,5,6,7,8,10,11'
failures=0

check() { # name type matrix gapO gapE targets queries
	local name=$1; shift
	"$WORK/ssw_native_check" "$@" 2>/dev/null | cut -f1-11 > "$WORK/$name.actual"
	zcat "$EXPECTED/$name.tsv.gz" | cut -f1-11 > "$WORK/$name.expected"
	local pairs
	pairs=$(wc -l < "$WORK/$name.expected")
	if cut -f$COMPARE "$WORK/$name.actual" | cmp -s - <(cut -f$COMPARE "$WORK/$name.expected"); then
		echo "PASS $name: $pairs pairs identical to the C reference"
	else
		echo "FAIL $name: differences (first ones below)"
		diff <(cut -f$COMPARE "$WORK/$name.actual") <(cut -f$COMPARE "$WORK/$name.expected") | head -10
		failures=$((failures + 1))
	fi
}

cat "$DATA/uniprot_human.fasta" "$DATA/synthetic_queries.fasta" > "$WORK/queries.fasta"
cat "$DATA/uniprot_homologs.fasta" "$DATA/synthetic_targets.fasta" > "$WORK/targets.fasta"

check protein_blosum50_o3_e1 protein blosum50 3 1 1 "$WORK/targets.fasta" "$WORK/queries.fasta"
check dna_m2_x2_o3_e1 dna 2:2 3 1 1 "$DATA/dna_1k.fa" "$DATA/dna_54mer_hap1_1.100.fa"

if [[ $failures -gt 0 ]]; then
	echo "$failures check(s) failed" >&2
	exit 1
fi
echo "all native checks passed"
