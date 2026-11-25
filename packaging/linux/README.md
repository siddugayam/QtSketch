# Linux build

The Linux package is a self-contained `.AppImage` (`x86_64` or `aarch64`) that runs on distributions with glibc 2.39 or newer (Ubuntu 24.04+, Debian 13+, Fedora 40+, Arch Linux).

## Option 1: One-command native setup & AppImage build

From the repository root on Ubuntu/Debian, Fedora/RHEL, Arch, or openSUSE:

```bash
./install_and_build.sh
```

- Installs required system build, X11/Wayland/OpenGL, and AppImage packages automatically.
- Uses the distribution's Qt 6 directly when it is 6.6 or newer, or automatically installs Qt 6.8.3 locally into `.qt/` via `aqtinstall` on older LTS releases (such as Ubuntu 22.04 and Ubuntu 24.04, which ship Qt 6.4 in `apt`).
- Runs `packaging/linux/fetch_indigo.sh` and `packaging/linux/build_imago.sh` (pass `--skip-imago` for a fast build without OCR).
- Builds `build_linux/QtSketch`, runs the 28 CTest suites headless (`QT_QPA_PLATFORM=offscreen`), and creates `dist/QtSketch-0.1-<arch>.AppImage` + `.sha256`.

## Option 2: Reproducible AppImage via Docker (any host)

```bash
packaging/linux/build_in_docker.sh amd64     # or arm64
```

Pass `--skip-imago` if you want to skip building Imago inside the container. The repository is mounted read-only and the finished `.AppImage` + `.sha256` are written to `dist/`.

To run the AppImage:

```bash
chmod +x QtSketch-0.1-x86_64.AppImage
./QtSketch-0.1-x86_64.AppImage
```

If FUSE 2 (`libfuse2` / `libfuse2t64`) is not installed, run with `--appimage-extract-and-run`.

## Option 3: Automated GitHub Actions / Releases

[`.github/workflows/release.yml`](../../.github/workflows/release.yml) builds macOS `arm64` `.dmg` and Linux `x86_64` + `aarch64` `.AppImage` packages in CI (caching Imago and Qt 6.8.3), uploads workflow artifacts on `workflow_dispatch`, and attaches all packages + SHA-256 files to a GitHub Release when a `v*` tag is pushed.
