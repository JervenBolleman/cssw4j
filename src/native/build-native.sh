#!/usr/bin/env bash
# Builds the cssw4j C library (libcssw4j.so and libcssw4j.h) with GraalVM Native Image.
#
#   src/native/build-native.sh <classes dir> <output dir>
#
# Needs a GraalVM JDK: set GRAALVM_HOME to its root (the directory containing bin/native-image). Without it the
# native-image found on PATH is used.
# Extra native-image options can be passed in NATIVE_IMAGE_ARGS, e.g. NATIVE_IMAGE_ARGS="-march=native".
# The classpath of the provided dependencies is read from <output dir>/classpath.txt, written by Maven.
set -euo pipefail

CLASSES="${1:?usage: build-native.sh <classes dir> <output dir>}"
OUT="${2:?usage: build-native.sh <classes dir> <output dir>}"
if [[ -n "${GRAALVM_HOME:-}" ]]; then
	NATIVE_IMAGE="$GRAALVM_HOME/bin/native-image"
else
	NATIVE_IMAGE="$(command -v native-image || true)"
fi

if [[ ! -x "$NATIVE_IMAGE" ]]; then
	echo "error: native-image not found. Set GRAALVM_HOME to a GraalVM JDK." >&2
	exit 1
fi

CP_FILE="${OUT%/}/classpath.txt"
DEPS=""
if [[ -f "$CP_FILE" ]]; then
	DEPS="$(cat "$CP_FILE")"
fi

mkdir -p "$OUT"
CLASSES="$(cd "$CLASSES" && pwd)"
OUT="$(cd "$OUT" && pwd)"
cd "$OUT"

# shellcheck disable=SC2086 # NATIVE_IMAGE_ARGS is a list of options
"$NATIVE_IMAGE" \
	--shared \
	--add-modules jdk.incubator.vector \
	-cp "$CLASSES${DEPS:+:$DEPS}" \
	-o libcssw4j \
	${NATIVE_IMAGE_ARGS:-} \
	swiss.sib.cssw4j.capi.CApi

echo "built: $OUT/libcssw4j.so"
echo "header: $OUT/libcssw4j.h"
