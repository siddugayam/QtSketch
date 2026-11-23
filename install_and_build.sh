#!/bin/bash
# One-command setup, build, test, and packaging script for QtSketch on macOS and Linux.
#
# Usage:
#   ./install_and_build.sh [options]
#
# Options:
#   --skip-imago     Skip building EPAM Imago from source (faster build; disables Insert Image OCR)
#   --skip-tests     Skip running the CTest test suite
#   --skip-package   Build the binary/app only; skip creating .dmg (macOS) or .AppImage (Linux)
#   --no-deps        Do not install system packages via brew / apt / dnf / pacman / zypper
#   --run            Launch QtSketch after the build completes
#   --jobs N         Override parallel build job count
#   -h, --help       Show this help message
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

SKIP_IMAGO=0
SKIP_TESTS=0
SKIP_PACKAGE=0
INSTALL_DEPS=1
RUN_AFTER=0
JOBS=""

while [ $# -gt 0 ]; do
    case "$1" in
        --skip-imago)   SKIP_IMAGO=1 ;;
        --skip-tests)   SKIP_TESTS=1 ;;
        --skip-package) SKIP_PACKAGE=1 ;;
        --no-deps)      INSTALL_DEPS=0 ;;
        --run)          RUN_AFTER=1 ;;
        --jobs)
            shift
            [ $# -gt 0 ] || { echo "--jobs requires an argument"; exit 1; }
            JOBS="$1"
            ;;
        -h|--help)
            sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1 (try --help)"
            exit 1
            ;;
    esac
    shift
done

run_privileged() {
    if [ "$(id -u)" -eq 0 ]; then
        "$@"
    elif command -v sudo >/dev/null 2>&1; then
        sudo "$@"
    else
        echo "Root privileges (or sudo) are required to install system packages: $*"
        exit 1
    fi
}

version_ge() {
    # Returns 0 if $1 >= $2 (dotted numeric versions).
    python3 - "$1" "$2" <<'PY'
import sys
a = [int(x) for x in sys.argv[1].split(".")[:3] if x.isdigit()]
b = [int(x) for x in sys.argv[2].split(".")[:3] if x.isdigit()]
sys.exit(0 if a >= b else 1)
PY
}

OS="$(uname -s)"

if [ "$OS" = "Darwin" ]; then
    echo "==> Detected macOS ($(uname -m))"
    if ! xcode-select -p >/dev/null 2>&1; then
        echo "Xcode Command Line Tools are required. Triggering installer..."
        xcode-select --install || true
        echo "Re-run ./install_and_build.sh once Xcode Command Line Tools finish installing."
        exit 1
    fi

    if [ "$INSTALL_DEPS" -eq 1 ]; then
        if ! command -v brew >/dev/null 2>&1; then
            echo "Homebrew is not installed. Install it from https://brew.sh and re-run, or pass --no-deps."
            exit 1
        fi
        for pkg in qt cmake ninja python; do
            if ! brew list "$pkg" >/dev/null 2>&1; then
                echo "==> Installing $pkg via Homebrew..."
                brew install "$pkg"
            fi
        done
    fi

    QT_PREFIX="${CMAKE_PREFIX_PATH:-$(brew --prefix qt)}"
    JOBS="${JOBS:-$(sysctl -n hw.logicalcpu)}"
    export JOBS

    echo "==> Fetching EPAM Indigo..."
    "$ROOT/packaging/macos/fetch_indigo.sh"

    WITH_IMAGO=ON
    if [ "$SKIP_IMAGO" -eq 1 ]; then
        WITH_IMAGO=OFF
        echo "==> Skipping Imago build (--skip-imago)"
    else
        echo "==> Building EPAM Imago (first run compiles OpenCV + Indigo static; cached thereafter)..."
        "$ROOT/packaging/macos/build_imago.sh"
    fi

    echo "==> Configuring and building QtSketch (build_mac)..."
    GEN_ARGS=()
    command -v ninja >/dev/null 2>&1 && GEN_ARGS=(-G Ninja)
    cmake -S "$ROOT" -B "$ROOT/build_mac" "${GEN_ARGS[@]}" \
        -DCMAKE_PREFIX_PATH="$QT_PREFIX" \
        -DCMAKE_BUILD_TYPE=Release \
        -DQTSKETCH_WITH_IMAGO="$WITH_IMAGO"
    cmake --build "$ROOT/build_mac" -j "$JOBS"

    if [ "$SKIP_TESTS" -eq 0 ]; then
        echo "==> Running CTest suite..."
        QT_QPA_PLATFORM=offscreen ctest --test-dir "$ROOT/build_mac" --output-on-failure -j "$JOBS"
    fi

    if [ "$SKIP_PACKAGE" -eq 0 ]; then
        echo "==> Packaging self-contained QtSketch.app and DMG..."
        QT_PREFIX="$QT_PREFIX" "$ROOT/packaging/macos/make_app.sh" "$ROOT/build_mac"
        "$ROOT/packaging/macos/make_dmg.sh"
        echo "==> Ready:"
        echo "    App: $ROOT/dist/QtSketch.app"
        echo "    DMG: $ROOT/dist/QtSketch-0.1-$(uname -m).dmg"
    fi

    if [ "$RUN_AFTER" -eq 1 ]; then
        if [ "$SKIP_PACKAGE" -eq 0 ] && [ -d "$ROOT/dist/QtSketch.app" ]; then
            open "$ROOT/dist/QtSketch.app"
        else
            open "$ROOT/build_mac/QtSketch.app"
        fi
    fi

elif [ "$OS" = "Linux" ]; then
    echo "==> Detected Linux ($(uname -m))"
    JOBS="${JOBS:-$(nproc)}"
    export JOBS

    if [ "$INSTALL_DEPS" -eq 1 ]; then
        if command -v apt-get >/dev/null 2>&1; then
            echo "==> Installing system dependencies via apt-get..."
            run_privileged apt-get update
            FUSE_PKG="libfuse2"
            apt-cache show libfuse2t64 >/dev/null 2>&1 && FUSE_PKG="libfuse2t64"
            run_privileged apt-get install -y \
                build-essential cmake ninja-build git curl ca-certificates file patchelf \
                desktop-file-utils squashfs-tools "$FUSE_PKG" \
                python3 python3-pip python3-venv \
                libgl1-mesa-dev libopengl0 libglx0 libegl1 \
                libxkbcommon0 libxkbcommon-x11-0 libfontconfig1 libfreetype6 \
                libdbus-1-3 libglib2.0-0 \
                libxcb-cursor0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 \
                libxcb-randr0 libxcb-render-util0 libxcb-shape0 libxcb-shm0 \
                libxcb-sync1 libxcb-xfixes0 libxcb-xinerama0 libxcb-xkb1 \
                libx11-xcb1 libxrender1 libxi6 libsm6 libice6 \
                libwayland-client0 libwayland-cursor0 libwayland-egl1
            # Install distro Qt 6 packages if available (used directly when >= 6.6).
            run_privileged apt-get install -y \
                qt6-base-dev qt6-declarative-dev qt6-svg-dev \
                qml6-module-qtquick-controls qml6-module-qtquick-layouts \
                qml6-module-qtquick-dialogs qml6-module-qtquick-shapes 2>/dev/null || true
        elif command -v dnf >/dev/null 2>&1; then
            echo "==> Installing system dependencies via dnf..."
            run_privileged dnf install -y \
                gcc gcc-c++ make cmake ninja-build git curl ca-certificates file patchelf \
                desktop-file-utils squashfs-tools fuse-libs python3 python3-pip \
                mesa-libGL-devel libxkbcommon-x11 fontconfig freetype dbus-libs glib2 \
                qt6-qtbase-devel qt6-qtdeclarative-devel qt6-qtsvg-devel qt6-qtwayland
        elif command -v pacman >/dev/null 2>&1; then
            echo "==> Installing system dependencies via pacman..."
            run_privileged pacman -Syu --needed --noconfirm \
                base-devel cmake ninja git curl ca-certificates file patchelf \
                desktop-file-utils squashfs-tools fuse2 python python-pip \
                qt6-base qt6-declarative qt6-svg qt6-wayland
        elif command -v zypper >/dev/null 2>&1; then
            echo "==> Installing system dependencies via zypper..."
            run_privileged zypper install -y \
                gcc gcc-c++ make cmake ninja git curl ca-certificates file patchelf \
                desktop-file-utils squashfs fuse python3 python3-pip \
                qt6-base-devel qt6-declarative-devel qt6-svg-devel
        else
            echo "Unrecognised Linux package manager; ensure CMake >= 3.24, C++20 compiler, Python 3, and Qt >= 6.6 are installed."
        fi
    fi

    # Ensure CMake >= 3.24 (e.g. Ubuntu 22.04 apt ships CMake 3.22, below CMakeLists.txt's 3.24 minimum).
    CMAKE_VER="$(cmake --version 2>/dev/null | awk 'NR==1{print $3}' || true)"
    if [ -z "$CMAKE_VER" ] || ! version_ge "$CMAKE_VER" "3.24.0"; then
        echo "==> System CMake is '${CMAKE_VER:-missing}' (< 3.24); installing CMake and Ninja into .aqt-venv..."
        if [ ! -d "$ROOT/.aqt-venv" ]; then
            python3 -m venv "$ROOT/.aqt-venv"
        fi
        "$ROOT/.aqt-venv/bin/pip" install -q --upgrade pip cmake ninja
        export PATH="$ROOT/.aqt-venv/bin:$PATH"
        echo "==> Using $(cmake --version | head -n 1)"
    fi

    # Detect Qt >= 6.6; if distro Qt is missing or older (e.g. Ubuntu 22.04 / 24.04 ships Qt 6.4),
    # install Qt 6.8.3 locally into .qt/ via aqtinstall.
    QT_PREFIX="${CMAKE_PREFIX_PATH:-}"
    QMAKE_BIN=""
    if [ -n "$QT_PREFIX" ] && [ -x "$QT_PREFIX/bin/qmake" ]; then
        QMAKE_BIN="$QT_PREFIX/bin/qmake"
    else
        QMAKE_BIN="$(command -v qmake6 || command -v qmake || true)"
    fi

    SYS_QT_VER=""
    if [ -n "$QMAKE_BIN" ]; then
        SYS_QT_VER="$("$QMAKE_BIN" -query QT_VERSION 2>/dev/null || true)"
    fi

    if [ -n "$SYS_QT_VER" ] && version_ge "$SYS_QT_VER" "6.6.0"; then
        if [ -z "$QT_PREFIX" ]; then
            QT_PREFIX="$("$QMAKE_BIN" -query QT_INSTALL_PREFIX 2>/dev/null || true)"
        fi
        echo "==> Using Qt $SYS_QT_VER ($QMAKE_BIN, prefix=${QT_PREFIX:-system})"
        export QMAKE="$QMAKE_BIN"
    else
        QT_VER="${QT_VERSION:-6.8.3}"
        echo "==> System Qt is '${SYS_QT_VER:-missing}' (< 6.6); installing Qt $QT_VER locally into .qt/ via aqtinstall..."
        if [ ! -d "$ROOT/.aqt-venv" ]; then
            python3 -m venv "$ROOT/.aqt-venv"
        fi
        "$ROOT/.aqt-venv/bin/pip" install -q --upgrade pip aqtinstall
        case "$(uname -m)" in
            aarch64|arm64) AQT_HOST="linux_arm64"; AQT_ARCH="linux_gcc_arm64" ;;
            *)             AQT_HOST="linux";       AQT_ARCH="linux_gcc_64" ;;
        esac
        if ! compgen -G "$ROOT/.qt/$QT_VER/gcc_*/bin/qmake" >/dev/null; then
            "$ROOT/.aqt-venv/bin/aqt" install-qt -O "$ROOT/.qt" "$AQT_HOST" desktop "$QT_VER" "$AQT_ARCH" -m qtshadertools
        fi
        QT_PREFIX=""
        for d in "$ROOT/.qt/$QT_VER"/gcc_*; do
            if [ -d "$d" ]; then
                QT_PREFIX="$d"
                break
            fi
        done
        [ -n "$QT_PREFIX" ] || { echo "Failed to locate installed Qt under $ROOT/.qt/$QT_VER"; exit 1; }
        export QMAKE="$QT_PREFIX/bin/qmake"
        export PATH="$QT_PREFIX/bin:$PATH"
        echo "==> Using local Qt $QT_VER at $QT_PREFIX"
    fi

    echo "==> Fetching EPAM Indigo..."
    "$ROOT/packaging/linux/fetch_indigo.sh"

    WITH_IMAGO=ON
    if [ "$SKIP_IMAGO" -eq 1 ]; then
        WITH_IMAGO=OFF
        echo "==> Skipping Imago build (--skip-imago)"
    else
        echo "==> Building EPAM Imago (first run compiles OpenCV + Indigo static; cached thereafter)..."
        "$ROOT/packaging/linux/build_imago.sh"
    fi

    echo "==> Configuring and building QtSketch (build_linux)..."
    GEN_ARGS=()
    command -v ninja >/dev/null 2>&1 && GEN_ARGS=(-G Ninja)
    CMAKE_QT_ARGS=()
    [ -n "$QT_PREFIX" ] && CMAKE_QT_ARGS=(-DCMAKE_PREFIX_PATH="$QT_PREFIX")
    cmake -S "$ROOT" -B "$ROOT/build_linux" "${GEN_ARGS[@]}" "${CMAKE_QT_ARGS[@]}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DQTSKETCH_WITH_IMAGO="$WITH_IMAGO"
    cmake --build "$ROOT/build_linux" -j "$JOBS"

    if [ "$SKIP_TESTS" -eq 0 ]; then
        echo "==> Running CTest suite..."
        QT_QPA_PLATFORM=offscreen ctest --test-dir "$ROOT/build_linux" --output-on-failure -j "$JOBS"
    fi

    if [ "$SKIP_PACKAGE" -eq 0 ]; then
        echo "==> Packaging self-contained AppImage..."
        "$ROOT/packaging/linux/make_appimage.sh" "$ROOT/build_linux"
        for appimg in "$ROOT"/dist/QtSketch-*.AppImage; do
            [ -f "$appimg" ] && echo "==> Ready: $appimg"
        done
    fi

    if [ "$RUN_AFTER" -eq 1 ]; then
        "$ROOT/build_linux/QtSketch" &
    fi
else
    echo "Unsupported OS '$OS'. On Windows, follow the PowerShell steps in README.md."
    exit 1
fi
