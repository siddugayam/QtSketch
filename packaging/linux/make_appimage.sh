#!/bin/bash
# Makes dist/QtSketch-<version>-<arch>.AppImage from a CMake build of the repo.
# Usage: packaging/linux/make_appimage.sh [build-dir]   (default: build)
# Needs qmake on PATH (or QMAKE set) for the Qt that the build used.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD="${1:-$ROOT/build}"
VERSION=0.1
ARCH="$(uname -m)"
APPDIR="$ROOT/dist/AppDir"
TOOLS="$ROOT/dist/linux-tools"

[ -x "$BUILD/QtSketch" ] || { echo "no $BUILD/QtSketch; build first"; exit 1; }
export QMAKE="${QMAKE:-$(command -v qmake6 || command -v qmake || true)}"
[ -n "$QMAKE" ] || { echo "qmake not found"; exit 1; }

mkdir -p "$TOOLS"
for t in linuxdeploy linuxdeploy-plugin-qt; do
    f="$TOOLS/$t-$ARCH.AppImage"
    [ -x "$f" ] || { curl -fsSL "https://github.com/linuxdeploy/$t/releases/download/continuous/$t-$ARCH.AppImage" -o "$f"; chmod +x "$f"; }
done
# Docker has no FUSE; the tools unpack themselves instead.
export APPIMAGE_EXTRACT_AND_RUN=1

rm -rf "$APPDIR"
mkdir -p "$APPDIR/usr/bin" "$APPDIR/usr/lib"
cp "$BUILD/QtSketch" "$APPDIR/usr/bin/"
cp "$ROOT"/indigo/lib/libindigo*.so "$APPDIR/usr/lib/"

# IndigoService looks for <exe dir>/../indigo and DocumentManager walks up to templates/.
mkdir -p "$APPDIR/usr/indigo/data/molecules/basic"
cp "$ROOT/assets/monomer_library.ket" "$APPDIR/usr/indigo/data/molecules/basic/"
cp -R "$ROOT/templates" "$APPDIR/usr/templates"

LIC="$APPDIR/usr/share/licenses/qtsketch"
mkdir -p "$LIC"
cp "$ROOT/LICENSE" "$ROOT/NOTICE" "$ROOT/THIRD-PARTY-NOTICES.md" "$LIC/"
cp "$ROOT/indigo/LICENCE-InChI" "$LIC/InChI-Licence-1.0.txt"
for l in LGPL-3.0-only GPL-3.0-only; do
    curl -fsSL --max-time 60 "https://raw.githubusercontent.com/spdx/license-list-data/main/text/$l.txt" -o "$LIC/$l.txt"
done

cp "$ROOT/packaging/linux/qtsketch.desktop" "$ROOT/dist/qtsketch.desktop"
cp "$ROOT/icons/app_icon.png" "$ROOT/dist/qtsketch.png"

export QML_SOURCES_PATHS="$ROOT"
export EXTRA_QT_MODULES="svg"
export LD_LIBRARY_PATH="$APPDIR/usr/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export LDAI_OUTPUT="$ROOT/dist/QtSketch-$VERSION-$ARCH.AppImage"
"$TOOLS/linuxdeploy-$ARCH.AppImage" --appdir "$APPDIR" \
    --desktop-file "$ROOT/dist/qtsketch.desktop" --icon-file "$ROOT/dist/qtsketch.png" \
    --plugin qt --output appimage
sha256sum "$LDAI_OUTPUT" | tee "$LDAI_OUTPUT.sha256"
