#!/bin/bash
# Builds EPAM's Imago (Apache-2.0, github.com/epam/Imago) from source and puts the
# resulting shared library and headers into <repo>/imago/.
#
# Unlike fetch_indigo.sh, this cannot just download a prebuilt package: EPAM
# publishes no Imago binaries or releases for any platform (confirmed: the
# GitHub repo has no Releases and no tags, and no PyPI wheel exists -- unlike
# epam.indigo, which fetch_indigo.sh downloads directly). Imago's own CMake
# build pulls its own vendored copies of Indigo and OpenCV as git submodules
# and builds all three from source, so this takes a while (OpenCV alone is a
# large build) and needs a C++ toolchain, cmake and git on PATH.
#
# The commit below is pinned to the last commit on the upstream repo (it has
# seen no activity since 2022) for a reproducible build.
set -euo pipefail
IMAGO_COMMIT=cd44e67d2692ed2aa82fa4db0df727adc6c77599
# Imago's own pinned third_party/opencv submodule commit at IMAGO_COMMIT, used
# only to fetch that exact OpenCV version's license text for attribution.
OPENCV_COMMIT=66f3c2673cd424a9c282fb0c8622a4ae43d3e81d
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DEST="$ROOT/imago"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

git clone https://github.com/epam/Imago.git "$TMP/src"
git -C "$TMP/src" checkout "$IMAGO_COMMIT"
git -C "$TMP/src" submodule update --init --recursive

cmake -S "$TMP/src" -B "$TMP/build" -DCMAKE_BUILD_TYPE=Release -DIMAGO_ENABLE_TESTS=OFF
cmake --build "$TMP/build" --config Release --target imago-c -- -j "$(nproc)"

mkdir -p "$DEST/lib" "$DEST/include"
cp "$TMP/build/api/c/libimago."* "$DEST/lib/"
cp "$ROOT/packaging/imago_c.h" "$DEST/include/imago_c.h"
cp "$TMP/src/LICENSE" "$DEST/LICENSE"
# OpenCV is statically linked into libimago -- its own license text travels with it.
curl -fsSL "https://raw.githubusercontent.com/opencv/opencv/$OPENCV_COMMIT/LICENSE" \
    -o "$DEST/LICENSE-OpenCV"
ls "$DEST/lib"
