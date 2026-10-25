#!/usr/bin/env bash
# Compiles the TornadoVM kernel module with the JDK the TornadoVM SDK runs on (JDK 21 for the "jdk21" SDKs, which are
# built with preview features; JDK 22 or newer for the "jdk22plus" SDKs).
#
#   TORNADO_JDK=/path/to/jdk TORNADOVM_HOME=/path/to/tornadovm tornado/build.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
: "${TORNADOVM_HOME:?set TORNADOVM_HOME to the TornadoVM SDK}"
JDK="${TORNADO_JDK:-${JAVA_HOME:?set TORNADO_JDK (or JAVA_HOME) to a JDK 21}}"
LIB="$TORNADOVM_HOME/share/java/tornado"
PREVIEW=()
if grep -qx 'tornado.jdk.preview=true' "$TORNADOVM_HOME/etc/tornado.jdk" 2>/dev/null; then
	PREVIEW=(--enable-preview)
fi
rm -rf "$HERE/target/classes" && mkdir -p "$HERE/target/classes"
"$JDK/bin/javac" "${PREVIEW[@]}" --release 21 \
	-cp "$LIB/tornado-api-7.3.1.jar:$LIB/tornado-annotation-7.3.1.jar" \
	-d "$HERE/target/classes" $(find "$HERE/src/main/java" -name '*.java')
echo "built $HERE/target/classes"
