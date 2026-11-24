#!/bin/bash
# Wraps dist/QtSketch.app in a compressed disk image with an Applications link. Run make_app.sh first.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
APP="$ROOT/dist/QtSketch.app"
OUT="$ROOT/dist/QtSketch-0.1-$(uname -m).dmg"
[ -d "$APP" ] || { echo "run make_app.sh first"; exit 1; }
codesign --verify --deep --strict "$APP"
S="$(mktemp -d)"; trap 'rm -rf "$S"' EXIT; cp -R "$APP" "$S/"; ln -s /Applications "$S/Applications"
rm -f "$OUT"
hdiutil create -volname QtSketch -srcfolder "$S" -ov -format ULMO "$OUT" >/dev/null 2>&1 \
    || hdiutil create -volname QtSketch -srcfolder "$S" -ov -format UDZO -imagekey zlib-level=9 "$OUT" >/dev/null
hdiutil verify "$OUT" >/dev/null
(cd "$ROOT/dist" && shasum -a 256 "$(basename "$OUT")" | tee "$(basename "$OUT").sha256")
