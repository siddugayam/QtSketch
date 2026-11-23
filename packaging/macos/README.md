# macOS build

Verified on Apple silicon (`arm64`) with Xcode Command Line Tools, CMake 4.x, and Homebrew Qt 6.

## One-command setup & DMG build

From the repository root:

```bash
./install_and_build.sh
```

Pass `--skip-imago` to skip building Imago from source if you do not need Insert Image OCR.

## Manual steps

```bash
brew install qt cmake ninja
packaging/macos/fetch_indigo.sh
packaging/macos/build_imago.sh          # optional; omit or pass -DQTSKETCH_WITH_IMAGO=OFF to skip
cmake -S . -B build_mac -G Ninja -DCMAKE_PREFIX_PATH="$(brew --prefix qt)" -DCMAKE_BUILD_TYPE=Release
cmake --build build_mac
QT_QPA_PLATFORM=offscreen ctest --test-dir build_mac --output-on-failure
packaging/macos/make_app.sh build_mac
packaging/macos/make_dmg.sh
```

- `fetch_indigo.sh` downloads the EPAM Indigo 1.45.0 wheel from PyPI using Python 3's standard library, verifies its SHA-256 digest, and places the libraries, C headers, and licences in `indigo/` (gitignored).
- `build_imago.sh` downloads pinned source archives for Imago, Indigo, and OpenCV, applies macOS SDK header fixes for bundled `zlib`/`libpng`, compiles only `opencv_core`, `opencv_imgproc`, `opencv_imgcodecs`, `indigo-static`, and `imago-c`, and installs `@rpath/libimago.dylib` into `imago/` (gitignored).
- `make_app.sh` bundles Qt, Indigo, Imago (if built), `monomer_library.ket`, `templates/`, and third-party licences inside `dist/QtSketch.app/Contents/`, rewrites any Homebrew load paths, and signs the bundle ad hoc.
- `make_dmg.sh` writes `dist/QtSketch-0.1-arm64.dmg` and its `.sha256` checksum.

## Limits

- The app is signed ad hoc, not notarised. Gatekeeper blocks it on first open; right-click the app and choose Open (or run `xattr -cr /Applications/QtSketch.app`).
- Some UI text specifies Windows font names first (`Segoe UI`, `Consolas`); macOS substitutes system fonts automatically.
