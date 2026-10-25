#!/usr/bin/env bash
# Runs VectorInlining with the vector helper inlined (default) and kept out of line; see README.md.
#   JAVA=/path/to/graalvm/bin/java ./run.sh
set -euo pipefail
cd "$(dirname "$0")"
JAVA="${JAVA:-java}"

echo "== default"
"$JAVA" --add-modules jdk.incubator.vector VectorInlining.java
echo "== shiftBytes not inlined"
"$JAVA" -XX:CompileCommand=quiet -XX:CompileCommand=dontinline,VectorInlining::shiftBytes \
	--add-modules jdk.incubator.vector VectorInlining.java
