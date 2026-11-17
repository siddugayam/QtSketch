#!/bin/bash
# Puts the EPAM Indigo 1.45.0 macOS arm64 libraries and C headers into <repo>/indigo/.
set -euo pipefail
VER=1.45.0
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DEST="$ROOT/indigo"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

python3 -m pip download "epam.indigo==$VER" --no-deps --only-binary=:all: \
    --platform macosx_11_0_arm64 -d "$TMP" >/dev/null
unzip -q "$TMP"/epam*indigo-*.whl -d "$TMP/whl"
mkdir -p "$DEST/lib"
cp "$TMP"/whl/indigo/lib/darwin-aarch64/libindigo*.dylib "$DEST/lib/"

RAW="https://raw.githubusercontent.com/epam/Indigo/indigo-$VER/api/c"
for h in indigo/indigo.h indigo-inchi/indigo-inchi.h indigo-renderer/indigo-renderer.h; do
    mkdir -p "$DEST/api/c/$(dirname "$h")"
    curl -fsSL "$RAW/$h" -o "$DEST/api/c/$h"
done
curl -fsSL "https://raw.githubusercontent.com/epam/Indigo/indigo-$VER/LICENSE" -o "$DEST/LICENSE"
# libindigo-inchi contains the InChI library, which has its own licence.
curl -fsSL "https://raw.githubusercontent.com/epam/Indigo/indigo-$VER/third_party/inchi/LICENCE" -o "$DEST/LICENCE-InChI"
ls "$DEST/lib"
