#!/bin/bash
# Builds, tests, and packages the Linux AppImage inside Docker from any host with Docker running.
# Usage: packaging/linux/build_in_docker.sh [amd64|arm64] [--skip-imago]
# Output: dist/QtSketch-0.1-<arch>.AppImage
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ARCH_ARG="amd64"
SKIP_IMAGO=0
for arg in "$@"; do
    case "$arg" in
        amd64|x86_64)  ARCH_ARG="amd64" ;;
        arm64|aarch64) ARCH_ARG="arm64" ;;
        --skip-imago)  SKIP_IMAGO=1 ;;
        *) echo "Usage: $0 [amd64|arm64] [--skip-imago]"; exit 1 ;;
    esac
done

PLAT="linux/$ARCH_ARG"
IMAGE="qtsketch-linux-build:$ARCH_ARG"

docker build --platform "$PLAT" -t "$IMAGE" "$ROOT/packaging/linux"
mkdir -p "$ROOT/dist"
# shellcheck disable=SC2016
docker run --rm --platform "$PLAT" -e SKIP_IMAGO="$SKIP_IMAGO" \
    -v "$ROOT:/src:ro" -v "$ROOT/dist:/out" "$IMAGE" bash -euo pipefail -c '
    mkdir /work && cd /src
    tar --exclude=./build --exclude=./build_* --exclude=./dist --exclude=./indigo \
        --exclude=./imago --exclude=./.qt --exclude=./.aqt-venv --exclude=./.git \
        -cf - . | tar -xf - -C /work
    cd /work
    packaging/linux/fetch_indigo.sh
    if [ "$SKIP_IMAGO" != "1" ]; then
        packaging/linux/build_imago.sh
    fi
    cmake -S . -B build_linux -G Ninja -DCMAKE_BUILD_TYPE=Release
    cmake --build build_linux
    QT_QPA_PLATFORM=offscreen ctest --test-dir build_linux --output-on-failure -j "$(nproc)"
    packaging/linux/make_appimage.sh build_linux
    cp dist/QtSketch-*.AppImage dist/QtSketch-*.AppImage.sha256 /out/
'
ls -lh "$ROOT"/dist/QtSketch-*.AppImage
