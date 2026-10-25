#!/usr/bin/env bash
# Runs the TornadoVM kernel on a device and compares it with the C reference results.
#
#   TORNADO_JDK=/path/to/jdk TORNADOVM_HOME=/path/to/tornadovm [LAYOUT=interleaved|contiguous] tornado/run-check.sh [jvm|tornado] [extra tornado options]
#
# "jvm" runs the same kernel as plain Java, without a device. "tornado" (default) needs a working OpenCL device:
# check `tornado --devices` first (the pocl package gives a CPU device).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
: "${TORNADOVM_HOME:?set TORNADOVM_HOME to the TornadoVM SDK}"
export JAVA_HOME="${TORNADO_JDK:-${JAVA_HOME:?set TORNADO_JDK (or JAVA_HOME) to a JDK 21}}"
MODE="${1:-tornado}"
shift || true
D="$ROOT/src/test/resources/dataset"
ARGS=("$D/blosum62.txt" "$D/uniprot_human.fasta,$D/synthetic_queries.fasta" "$D/uniprot_homologs.fasta,$D/synthetic_targets.fasta"
	11 1 "$ROOT/src/test/resources/expected/protein_blosum62_o11_e1.tsv.gz" "${LAYOUT:-interleaved}")
if [[ "$MODE" == "jvm" ]]; then
	LIB="$TORNADOVM_HOME/share/java/tornado"
	PREVIEW=()
	if grep -qx 'tornado.jdk.preview=true' "$TORNADOVM_HOME/etc/tornado.jdk" 2>/dev/null; then
		PREVIEW=(--enable-preview)
	fi
	exec "$JAVA_HOME/bin/java" "${PREVIEW[@]}" \
		-cp "$HERE/target/classes:$LIB/tornado-api-7.3.1.jar:$LIB/tornado-annotation-7.3.1.jar" \
		swiss.sib.cssw4j.tornado.Check jvm "${ARGS[@]}"
fi
exec "$TORNADOVM_HOME/bin/tornado" "$@" --classpath "$HERE/target/classes" \
	swiss.sib.cssw4j.tornado.Check tornado "${ARGS[@]}"
