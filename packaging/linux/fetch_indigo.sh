#!/bin/bash
# Puts the EPAM Indigo 1.45.0 Linux libraries (x86_64 or aarch64) and C headers into <repo>/indigo/.
set -euo pipefail
VER=1.45.0
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DEST="$ROOT/indigo"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

case "$(uname -m)" in
    x86_64)        PLAT=manylinux1_x86_64;     DIR=linux-x86_64 ;;
    aarch64|arm64) PLAT=manylinux2014_aarch64; DIR=linux-aarch64 ;;
    *) echo "unsupported CPU: $(uname -m)"; exit 1 ;;
esac

python3 -m pip download "epam.indigo==$VER" --no-deps --only-binary=:all: \
    --platform "$PLAT" -d "$TMP" >/dev/null
python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" "$TMP"/epam*indigo-*.whl "$TMP/whl"
mkdir -p "$DEST/lib"
cp "$TMP"/whl/indigo/lib/$DIR/libindigo*.so "$DEST/lib/"

RAW="https://raw.githubusercontent.com/epam/Indigo/indigo-$VER/api/c"
for h in indigo/indigo.h indigo-inchi/indigo-inchi.h indigo-renderer/indigo-renderer.h; do
    mkdir -p "$DEST/api/c/$(dirname "$h")"
    curl -fsSL "$RAW/$h" -o "$DEST/api/c/$h"
done
curl -fsSL "https://raw.githubusercontent.com/epam/Indigo/indigo-$VER/LICENSE" -o "$DEST/LICENSE"
# libindigo-inchi contains the InChI library, which has its own licence.
curl -fsSL "https://raw.githubusercontent.com/epam/Indigo/indigo-$VER/third_party/inchi/LICENCE" -o "$DEST/LICENCE-InChI"
ls "$DEST/lib"
