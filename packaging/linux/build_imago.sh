#!/bin/bash
# Builds EPAM Imago (Apache-2.0, github.com/epam/Imago) from pinned source archives into <repo>/imago/.
set -euo pipefail

IMAGO_COMMIT=cd44e67d2692ed2aa82fa4db0df727adc6c77599
INDIGO_SUBMODULE_COMMIT=fdb66e24526b698205e13d96a652cd4847c5e43b
OPENCV_SUBMODULE_COMMIT=66f3c2673cd424a9c282fb0c8622a4ae43d3e81d

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DEST="$ROOT/imago"

if [ "${1:-}" != "--force" ] \
    && [ -f "$DEST/lib/libimago.so" ] \
    && [ -f "$DEST/include/imago_c.h" ] \
    && [ -f "$DEST/LICENSE" ] \
    && [ -f "$DEST/LICENSE-OpenCV" ]; then
    echo "Imago already built in $DEST (pass --force to rebuild)"
    exit 0
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

mkdir -p "$TMP/src/third_party/indigo" "$TMP/src/third_party/opencv"

echo "Downloading pinned Imago, Indigo, and OpenCV source archives..."
curl -fsSL "https://codeload.github.com/epam/Imago/tar.gz/$IMAGO_COMMIT" \
    | tar -xz -C "$TMP/src" --strip-components=1
curl -fsSL "https://codeload.github.com/epam/Indigo/tar.gz/$INDIGO_SUBMODULE_COMMIT" \
    | tar -xz -C "$TMP/src/third_party/indigo" --strip-components=1
curl -fsSL "https://codeload.github.com/opencv/opencv/tar.gz/$OPENCV_SUBMODULE_COMMIT" \
    | tar -xz -C "$TMP/src/third_party/opencv" --strip-components=1

GEN_ARGS=()
if command -v ninja >/dev/null 2>&1; then
    GEN_ARGS=(-G Ninja)
fi

# Pre-include <cstdint> so 2021 OpenCV/Indigo sources compile cleanly on GCC 13/14/15+.
cmake -S "$TMP/src" -B "$TMP/build" "${GEN_ARGS[@]}" \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-include cstdint" \
    -DIMAGO_ENABLE_TESTS=OFF \
    -DBUILD_LIST=core,imgproc,imgcodecs \
    -DBUILD_TESTS=OFF \
    -DBUILD_PERF_TESTS=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_opencv_apps=OFF \
    -DBUILD_JAVA=OFF \
    -DBUILD_opencv_python2=OFF \
    -DBUILD_opencv_python3=OFF \
    -DWITH_OPENCL=OFF \
    -DWITH_IPP=OFF \
    -DWITH_ITT=OFF \
    -DWITH_PROTOBUF=OFF \
    -DWITH_WEBP=OFF \
    -DWITH_TIFF=OFF \
    -DWITH_JASPER=OFF \
    -DWITH_OPENJPEG=OFF \
    -DWITH_OPENEXR=OFF \
    -DWITH_ADE=OFF \
    -DWITH_1394=OFF \
    -DWITH_GSTREAMER=OFF \
    -DWITH_FFMPEG=OFF \
    -DWITH_V4L=OFF \
    -DPNG_ARM_NEON=off

JOBS="${JOBS:-$(nproc)}"
cmake --build "$TMP/build" --config Release --target imago-c -j "$JOBS"

mkdir -p "$DEST/lib" "$DEST/include"
if compgen -G "$TMP/build/bin/libimago.*" >/dev/null; then
    cp "$TMP/build/bin/libimago."* "$DEST/lib/"
else
    cp "$TMP/build/api/c/libimago."* "$DEST/lib/"
fi

cp "$ROOT/packaging/imago_c.h" "$DEST/include/imago_c.h"
cp "$TMP/src/LICENSE" "$DEST/LICENSE"
cp "$TMP/src/third_party/opencv/LICENSE" "$DEST/LICENSE-OpenCV"
ls -lh "$DEST/lib"
