#!/usr/bin/env bash
# Runs cssw4j's benchmarks for one or two builds of cssw4j, interleaved, and prints the results as Markdown.
#
# Usage: run-benchmarks.sh OUT CODE_DIR [REFERENCE_DIR]
#
# Each directory is a cssw4j checkout built with `mvn test-compile` and
# `mvn dependency:build-classpath -Dmdep.outputFile=target/cp.txt`. For each round, each variant and each tree (the
# trees innermost, so that code and reference see the same machine state) it runs
#   1. the dataset benchmark, swiss.sib.cssw4j.Benchmark: the protein test set all against all, scores only (flag 0)
#      and with CIGAR (flag 1), in GCUPS. This is the library's real aligner end to end.
#   2. the JMH benchmarks matching BENCHMARK.
# Every JVM gets BASE_OPTS and then the variant's own options. A variant whose options the JVM rejects (such as
# -XX:UseAVX=2 on AArch64) is skipped on this machine.
#
# Settings (environment):
#   VARIANTS      label=JVM options, separated by ;          (default: "default=")
#   BENCH_ROUNDS  rounds, the variants and trees interleaved  (default: 1)
#   BENCHMARK     JMH benchmarks, a regular expression        (default: .)
#   BENCH_ARGS    further JMH arguments, such as -f 1 -wi 1 -i 3
#   DATASET_ARGS  arguments of swiss.sib.cssw4j.Benchmark: rounds, gap open, gap extend (default: 20 10 2)
#   JAVA          the java launcher (default: java on the PATH)
#
# Writes to OUT: jmh-<tree>-<variant>-<round>.csv and .log, dataset-<tree>-<variant>-<round>.log, and dataset.tsv
# (one line per tree, variant, round and flag; read by report.py). Tree 0 is the code, tree 1 the reference.
set -euo pipefail

out=$(mkdir -p "$1" && cd "$1" && pwd)
shift
scripts=$(cd "$(dirname "$0")" && pwd)
JAVA=${JAVA:-java}
VARIANTS=${VARIANTS:-default=}
BENCH_ROUNDS=${BENCH_ROUNDS:-1}
BENCHMARK=${BENCHMARK:-.}
BENCH_ARGS=${BENCH_ARGS:-}
DATASET_ARGS=${DATASET_ARGS:-20 10 2}
BASE_OPTS="--add-modules jdk.incubator.vector"

# The trees: "label short-sha=dir". The code is labelled with its branch, or "code" for a detached checkout.
trees=()
for dir in "$@"; do
	dir=$(cd "$dir" && pwd)
	sha=$(git -C "$dir" rev-parse --short HEAD)
	if [ ${#trees[@]} -eq 0 ]; then
		label=$(git -C "$dir" rev-parse --abbrev-ref HEAD | sed 's/^HEAD$/code/')
	else
		label=reference
	fi
	trees+=("$label $sha=$dir")
done

# The variants this JVM accepts.
IFS=';' read -ra all <<<"$VARIANTS"
variants=()
skipped=()
for v in "${all[@]}"; do
	[ -n "$v" ] || continue
	[[ $v == *=* ]] || v="$v="
	# shellcheck disable=SC2086
	if "$JAVA" $BASE_OPTS ${v#*=} -version >/dev/null 2>&1; then
		variants+=("$v")
	else
		skipped+=("$v")
	fi
done

printf 'tree\tlabel\tvariant\tround\tlanes\tflag\tseconds\tgcups\tchecksum\n' >"$out/dataset.tsv"
for r in $(seq 1 "$BENCH_ROUNDS"); do
	for v in "${variants[@]}"; do
		name=${v%%=*}
		opts=${v#*=}
		for k in "${!trees[@]}"; do
			dir=${trees[$k]#*=}
			cp="$dir/target/test-classes:$dir/target/classes:$(cat "$dir/target/cp.txt")"
			log="$out/dataset-$k-$name-$r.log"
			# shellcheck disable=SC2086
			"$JAVA" $BASE_OPTS $opts -cp "$cp" swiss.sib.cssw4j.Benchmark $DATASET_ARGS >"$log" 2>&1 || true
			# Each flag is printed twice, warm-up first; the last line of each flag is the measurement.
			awk -v k="$k" -v label="${trees[$k]%%=*}" -v name="$name" -v r="$r" '
				/byte lanes/ { lanes = $(NF - 2) }
				/flag=[0-9]+ .* GCUPS/ {
					match($0, /flag=[0-9]+/); flag = substr($0, RSTART + 5, RLENGTH - 5)
					match($0, /: *[0-9.]+ s/); s = substr($0, RSTART, RLENGTH); gsub(/[: s]/, "", s)
					match($0, /[0-9.]+ GCUPS/); g = substr($0, RSTART, RLENGTH); sub(/ GCUPS/, "", g)
					match($0, /checksum -?[0-9]+/); c = substr($0, RSTART + 9, RLENGTH - 9)
					last[flag] = k "\t" label "\t" name "\t" r "\t" lanes "\t" flag "\t" s "\t" g "\t" c
					if (!(flag in seen)) { seen[flag] = 1; order[++n] = flag }
				}
				END { for (i = 1; i <= n; i++) print last[order[i]] }' "$log" >>"$out/dataset.tsv"
			# shellcheck disable=SC2086
			(cd "$dir" &&
				"$JAVA" $BASE_OPTS -cp "$cp" org.openjdk.jmh.Main "$BENCHMARK" $BENCH_ARGS \
					-jvmArgsAppend "$BASE_OPTS $opts" \
					-rf csv -rff "$out/jmh-$k-$name-$r.csv" >"$out/jmh-$k-$name-$r.log" 2>&1) || true
		done
	done
done

# The Markdown: the variants, the dataset table and the JMH table, one column per tree, variant and round.
echo
echo "### Variants"
echo
echo "Every JVM: \`$BASE_OPTS\`, then the variant's options."
echo
for v in "${variants[@]}"; do
	lanes=$(awk -F'\t' -v n="${v%%=*}" '$3 == n && $1 == 0 { print $5; exit }' "$out/dataset.tsv")
	opts=${v#*=}
	echo "- ${v%%=*}: ${opts:+\`$opts\`}${opts:-no further options}${lanes:+, $lanes byte lanes}"
done
for v in "${skipped[@]}"; do
	echo "- ${v%%=*}: \`${v#*=}\`, skipped: this JVM does not accept the options"
done
echo
echo "### Dataset: protein set all against all ($DATASET_ARGS: rounds, gap open, gap extend), GCUPS, one thread"
echo
python3 "$scripts/report.py" dataset "$out/dataset.tsv"
runs=()
for r in $(seq 1 "$BENCH_ROUNDS"); do
	for v in "${variants[@]}"; do
		for k in "${!trees[@]}"; do
			label="${trees[$k]%%=*} ${v%%=*}"
			[ "$BENCH_ROUNDS" -gt 1 ] && label+=" $r"
			runs+=("$label=$out/jmh-$k-${v%%=*}-$r.csv")
		done
	done
done
echo
echo "### JMH: $BENCHMARK $BENCH_ARGS (score ± 99.9% error)"
echo
python3 "$scripts/jmh-table.py" "${runs[@]}"
