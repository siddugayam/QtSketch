#!/bin/bash
# Builds, tests and packages the Linux AppImage inside Docker, from any host that has Docker.
# Usage: packaging/linux/build_in_docker.sh [amd64|arm64]   (default: amd64)
# Output: dist/QtSketch-0.1-<arch>.AppImage. The repo itself is only read; work happens on a copy.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PLAT="linux/${1:-amd64}"
IMAGE="qtsketch-linux-build:${1:-amd64}"

docker build --platform "$PLAT" -t "$IMAGE" "$ROOT/packaging/linux"
mkdir -p "$ROOT/dist"
docker run --rm --platform "$PLAT" -v "$ROOT:/src:ro" -v "$ROOT/dist:/out" "$IMAGE" bash -euo pipefail -c '
    mkdir /work && cd /src
    tar --exclude=./build --exclude=./build_* --exclude=./dist --exclude=./indigo --exclude=./.git -cf - . | tar -xf - -C /work
    cd /work
    packaging/linux/fetch_indigo.sh
    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
    cmake --build build
    QT_QPA_PLATFORM=offscreen ctest --test-dir build --output-on-failure
    packaging/linux/make_appimage.sh build
    cp dist/QtSketch-*.AppImage dist/QtSketch-*.AppImage.sha256 /out/
'
ls -l "$ROOT"/dist/QtSketch-*.AppImage
