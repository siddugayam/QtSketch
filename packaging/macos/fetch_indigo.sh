#!/bin/bash
# Downloads EPAM Indigo 1.45.0 macOS shared libraries, C headers, and license files into <repo>/indigo/.
# Uses only Python 3's standard library and curl (no pip or unzip required).
set -euo pipefail

VER="${INDIGO_VERSION:-1.45.0}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DEST="$ROOT/indigo"

if [ "${1:-}" != "--force" ] \
    && [ -f "$DEST/lib/libindigo.dylib" ] \
    && [ -f "$DEST/lib/libindigo-inchi.dylib" ] \
    && [ -f "$DEST/lib/libindigo-renderer.dylib" ] \
    && [ -f "$DEST/api/c/indigo/indigo.h" ] \
    && [ -f "$DEST/LICENSE" ] \
    && [ -f "$DEST/LICENCE-InChI" ]; then
    echo "Indigo $VER already present in $DEST (pass --force to re-download)"
    exit 0
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

ARCH="$(uname -m)"
mkdir -p "$DEST/lib"

python3 - "$VER" "$ARCH" "$TMP" "$DEST/lib" <<'PY'
import hashlib
import json
import os
import shutil
import sys
import urllib.request
import zipfile

ver, arch, tmp_dir, out_lib = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
url = f"https://pypi.org/pypi/epam.indigo/{ver}/json"
req = urllib.request.Request(url, headers={"User-Agent": "QtSketch-fetch-indigo/1.0"})
with urllib.request.urlopen(req, timeout=60) as resp:  # nosec B310
    meta = json.load(resp)

tags = ("arm64", "universal2") if arch == "arm64" else ("x86_64", "universal2")
chosen = None
for item in meta.get("urls", []):
    fn = item.get("filename", "")
    if fn.endswith(".whl") and "macosx" in fn and any(t in fn for t in tags):
        chosen = item
        break

if not chosen:
    sys.exit(f"No macOS {arch} wheel found for epam.indigo=={ver}")
if not chosen["url"].startswith("https://"):
    sys.exit(f"Refusing non-HTTPS wheel URL: {chosen['url']}")

whl_path = os.path.join(tmp_dir, chosen["filename"])
req_w = urllib.request.Request(chosen["url"], headers={"User-Agent": "QtSketch-fetch-indigo/1.0"})
with urllib.request.urlopen(req_w, timeout=120) as resp, open(whl_path, "wb") as out:  # nosec B310
    shutil.copyfileobj(resp, out)

expected_sha = chosen.get("digests", {}).get("sha256")
if expected_sha:
    with open(whl_path, "rb") as fh:
        actual_sha = hashlib.sha256(fh.read()).hexdigest()
    if actual_sha != expected_sha:
        sys.exit(f"SHA-256 mismatch for {chosen['filename']}: {actual_sha} != {expected_sha}")

copied = 0
with zipfile.ZipFile(whl_path) as zf:
    for name in zf.namelist():
        base = os.path.basename(name)
        if "/lib/" in name and base.startswith("libindigo") and base.endswith(".dylib"):
            with zf.open(name) as src, open(os.path.join(out_lib, base), "wb") as dst:
                shutil.copyfileobj(src, dst)
            copied += 1

if copied == 0:
    sys.exit(f"No libindigo*.dylib found inside {chosen['filename']}")
PY

RAW="https://raw.githubusercontent.com/epam/Indigo/indigo-$VER/api/c"
for h in indigo/indigo.h indigo-inchi/indigo-inchi.h indigo-renderer/indigo-renderer.h; do
    mkdir -p "$DEST/api/c/$(dirname "$h")"
    curl -fsSL "$RAW/$h" -o "$DEST/api/c/$h"
done
curl -fsSL "https://raw.githubusercontent.com/epam/Indigo/indigo-$VER/LICENSE" -o "$DEST/LICENSE"
curl -fsSL "https://raw.githubusercontent.com/epam/Indigo/indigo-$VER/third_party/inchi/LICENCE" -o "$DEST/LICENCE-InChI"
ls "$DEST/lib"
