# Linux build

The Linux package is an AppImage: one file that runs on most current distributions without
installing anything. It is built inside Docker on Ubuntu 24.04 with Qt 6.8 LTS from qt.io, so it
needs glibc 2.39 or newer (Ubuntu 24.04, Debian 13, Fedora 40 and later).

## AppImage with Docker (any host)

```bash
packaging/linux/build_in_docker.sh amd64     # or arm64
```

This builds the image in `packaging/linux/Dockerfile`, copies the repo into the container, fetches
Indigo, builds, runs the tests with the offscreen Qt platform, and writes
`dist/QtSketch-0.1-x86_64.AppImage` (or `-aarch64`) with its SHA-256. The repo itself is mounted
read-only.

To run it:

```bash
chmod +x QtSketch-0.1-x86_64.AppImage
./QtSketch-0.1-x86_64.AppImage
```

If the system has no FUSE 2 (`libfuse2`), run it with `--appimage-extract-and-run`.

## Native build (for development)

On a distribution with Qt 6.6 or newer. These package names have not been tested; the Docker
route above is the tested one.

```bash
# Debian 13 / Ubuntu 25.04 or later
sudo apt install build-essential cmake ninja-build python3-pip \
    qt6-base-dev qt6-declarative-dev qt6-svg-dev qml6-module-qtquick-controls \
    qml6-module-qtquick-layouts qml6-module-qtquick-dialogs qml6-module-qtquick-shapes
packaging/linux/fetch_indigo.sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
QT_QPA_PLATFORM=offscreen ctest --test-dir build
./build/QtSketch
```

- `fetch_indigo.sh` puts the Indigo 1.45.0 libraries, headers and licences in `indigo/`
  (gitignored). The libraries come from the `epam.indigo` wheel on PyPI for x86_64 or aarch64.
- `make_appimage.sh <build-dir>` bundles a build into an AppImage with linuxdeploy and its Qt
  plugin, which it downloads into `dist/linux-tools/`.

## Limits

- Image recognition (Imago) is not available. `ImagoServiceStub.cpp` returns an error instead.
- The AppImage is not signed.
