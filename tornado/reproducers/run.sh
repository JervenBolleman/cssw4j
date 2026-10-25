#!/usr/bin/env bash
# Compiles and runs the reproducer with a TornadoVM SDK.
#
#   TORNADOVM_HOME=/path/to/tornadovm JAVA_HOME=/path/to/jdk ./run.sh [device]   # device as driver:index, default 0:0
#
# With the CUDA SDK also set CUDA_HOME and LD_LIBRARY_PATH to the CUDA toolkit.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
: "${TORNADOVM_HOME:?set TORNADOVM_HOME}"
: "${JAVA_HOME:?set JAVA_HOME}"
LIB="$TORNADOVM_HOME/share/java/tornado"
API="$(ls "$LIB"/tornado-api-*.jar | head -1)"
ANN="$(ls "$LIB"/tornado-annotation-*.jar | head -1)"
mkdir -p "$HERE/classes"
"$JAVA_HOME/bin/javac" -cp "$API:$ANN" -d "$HERE/classes" "$HERE/ByteToInt.java"
DEVICE="${1:-0:0}"
exec "$TORNADOVM_HOME/bin/tornado" "--jvm=-Ds.t.device=$DEVICE" --classpath "$HERE/classes" ByteToInt
