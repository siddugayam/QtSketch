#!/bin/bash
# Builds dist/QtSketch-<version>-<arch>.AppImage from a CMake build of the repo.
# Usage: packaging/linux/make_appimage.sh [build-dir]   (default: build_linux, fallback: build)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
if [ $# -ge 1 ]; then
    BUILD="$1"
elif [ -x "$ROOT/build_linux/QtSketch" ]; then
    BUILD="$ROOT/build_linux"
else
    BUILD="$ROOT/build"
fi
case "$BUILD" in
    /*) ;;
    *) BUILD="$ROOT/$BUILD" ;;
esac

VERSION="${QTSKETCH_VERSION:-0.1}"
RAW_ARCH="$(uname -m)"
case "$RAW_ARCH" in
    x86_64|amd64)  ARCH="x86_64" ;;
    aarch64|arm64) ARCH="aarch64" ;;
    *) echo "Unsupported AppImage architecture: $RAW_ARCH"; exit 1 ;;
esac

APPDIR="$ROOT/dist/AppDir"
TOOLS="$ROOT/dist/linux-tools"

[ -x "$BUILD/QtSketch" ] || { echo "No $BUILD/QtSketch binary found; build first."; exit 1; }
export QMAKE="${QMAKE:-$(command -v qmake6 || command -v qmake || true)}"
[ -n "$QMAKE" ] || { echo "qmake6/qmake not found; put Qt's bin directory on PATH or set QMAKE."; exit 1; }

mkdir -p "$TOOLS" "$ROOT/dist"
for t in linuxdeploy linuxdeploy-plugin-qt; do
    f="$TOOLS/$t-$ARCH.AppImage"
    if [ ! -x "$f" ]; then
        curl -fsSL "https://github.com/linuxdeploy/$t/releases/download/continuous/$t-$ARCH.AppImage" -o "$f"
        chmod +x "$f"
    fi
done

# Extract-and-run avoids requiring FUSE inside containers or minimal hosts.
export APPIMAGE_EXTRACT_AND_RUN=1

rm -rf "$APPDIR"
mkdir -p "$APPDIR/usr/bin" "$APPDIR/usr/lib" "$APPDIR/usr/share/qtsketch"
cp "$BUILD/QtSketch" "$APPDIR/usr/bin/"
cp "$ROOT"/indigo/lib/libindigo*.so* "$APPDIR/usr/lib/"
if compgen -G "$ROOT/imago/lib/libimago*.so*" >/dev/null; then
    cp "$ROOT"/imago/lib/libimago*.so* "$APPDIR/usr/lib/"
fi

if command -v patchelf >/dev/null 2>&1; then
    # shellcheck disable=SC2016
    patchelf --set-rpath '$ORIGIN/../lib' "$APPDIR/usr/bin/QtSketch" || true
    for so in "$APPDIR/usr/lib"/libindigo*.so* "$APPDIR/usr/lib"/libimago*.so*; do
        # shellcheck disable=SC2016
        [ -f "$so" ] && patchelf --set-rpath '$ORIGIN' "$so" 2>/dev/null || true
    done
fi

cp "$ROOT/assets/monomer_library.ket" "$APPDIR/usr/share/qtsketch/monomer_library.ket"
cp -R "$ROOT/templates" "$APPDIR/usr/share/qtsketch/templates"
mkdir -p "$APPDIR/usr/indigo/data/molecules/basic"
cp "$ROOT/assets/monomer_library.ket" "$APPDIR/usr/indigo/data/molecules/basic/"
cp -R "$ROOT/templates" "$APPDIR/usr/templates"

LIC="$APPDIR/usr/share/licenses/qtsketch"
mkdir -p "$LIC"
for f in LICENSE NOTICE THIRD-PARTY-NOTICES.md; do
    [ -f "$ROOT/$f" ] && cp "$ROOT/$f" "$LIC/"
done
[ -f "$ROOT/indigo/LICENSE" ] && cp "$ROOT/indigo/LICENSE" "$LIC/Indigo-Apache-2.0.txt"
[ -f "$ROOT/indigo/LICENCE-InChI" ] && cp "$ROOT/indigo/LICENCE-InChI" "$LIC/InChI-Licence-1.0.txt"
[ -f "$ROOT/imago/LICENSE" ] && cp "$ROOT/imago/LICENSE" "$LIC/Imago-Apache-2.0.txt"
[ -f "$ROOT/imago/LICENSE-OpenCV" ] && cp "$ROOT/imago/LICENSE-OpenCV" "$LIC/OpenCV-Apache-2.0.txt"
for l in LGPL-3.0-only GPL-3.0-only; do
    [ -f "$LIC/$l.txt" ] || curl -fsSL --max-time 20 \
        "https://raw.githubusercontent.com/spdx/license-list-data/main/text/$l.txt" \
        -o "$LIC/$l.txt" 2>/dev/null || true
done

cp "$ROOT/packaging/linux/qtsketch.desktop" "$ROOT/dist/qtsketch.desktop"
cp "$ROOT/icons/app_icon.png" "$ROOT/dist/qtsketch.png"

export QML_SOURCES_PATHS="$ROOT"
export EXTRA_QT_MODULES="svg"
export LD_LIBRARY_PATH="$APPDIR/usr/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export LDAI_OUTPUT="$ROOT/dist/QtSketch-$VERSION-$ARCH.AppImage"
rm -f "$LDAI_OUTPUT" "$LDAI_OUTPUT.sha256"

"$TOOLS/linuxdeploy-$ARCH.AppImage" --appdir "$APPDIR" \
    --desktop-file "$ROOT/dist/qtsketch.desktop" \
    --icon-file "$ROOT/dist/qtsketch.png" \
    --plugin qt \
    --output appimage

(cd "$(dirname "$LDAI_OUTPUT")" && sha256sum "$(basename "$LDAI_OUTPUT")" | tee "$(basename "$LDAI_OUTPUT").sha256")
echo "Built $LDAI_OUTPUT"
