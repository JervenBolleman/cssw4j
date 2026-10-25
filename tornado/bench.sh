#!/usr/bin/env bash
# Runs the JMH benchmark of the kernel (src/bench) through a TornadoVM SDK, on the device named by cssw4j.device.
#
#   TORNADOVM_HOME=... TORNADO_JDK=... [DEVICE=0:0] [BENCH_FILTER=regex] [CUDA_HOME=... LD_LIBRARY_PATH=...] tornado/bench.sh
#
# Needs the JMH 1.37 jars in the local Maven repository (~/.m2), as fetched for the main build.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
: "${TORNADOVM_HOME:?set TORNADOVM_HOME to the TornadoVM SDK}"
JDK="${TORNADO_JDK:?set TORNADO_JDK to a JDK 22 or newer}"
M2="${MAVEN_REPO:-$HOME/.m2/repository}"
JMH="$M2/org/openjdk/jmh/jmh-core/1.37/jmh-core-1.37.jar"
ANN="$M2/org/openjdk/jmh/jmh-generator-annprocess/1.37/jmh-generator-annprocess-1.37.jar"
OPTS="$M2/net/sf/jopt-simple/jopt-simple/5.0.4/jopt-simple-5.0.4.jar"
MATH="$M2/org/apache/commons/commons-math3/3.2/commons-math3-3.2.jar"
for f in "$JMH" "$ANN" "$OPTS" "$MATH"; do
	[[ -f "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
LIB="$TORNADOVM_HOME/share/java/tornado"
TORNADO_CP="$LIB/tornado-api-7.3.1.jar:$LIB/tornado-annotation-7.3.1.jar"
JMH_CP="$JMH:$OPTS:$MATH"

# Compile the benchmark with the JMH annotation processor (generates the benchmark glue code).
rm -rf "$HERE/target/bench" && mkdir -p "$HERE/target/bench"
"$JDK/bin/javac" --release 21 -cp "$HERE/target/classes:$TORNADO_CP:$JMH_CP" \
	-processorpath "$ANN:$JMH_CP" -d "$HERE/target/bench" \
	$(find "$HERE/src/bench/java" -name '*.java')

# DEVICE selects the TornadoVM device (driver:index, from `tornado --devices`); the launcher adds the TornadoVM flags.
cd "$HERE"
export JAVA_HOME="$JDK" # the launcher runs the JVM from JAVA_HOME
DEVICE_OPT=()
if [[ -n "${DEVICE:-}" ]]; then
	DEVICE_OPT=("--jvm=-Dcssw4j.device=$DEVICE")
fi
exec "$TORNADOVM_HOME/bin/tornado" "${DEVICE_OPT[@]}" \
	--classpath "$HERE/target/classes:$HERE/target/bench:$JMH_CP" \
	swiss.sib.cssw4j.tornado.KernelBenchRunner "${BENCH_FILTER:-KernelBench}"
