# macOS build

Tested on an Apple silicon Mac with macOS 26, Xcode command line tools, Homebrew Qt 6.12 and CMake.

## Steps

```bash
brew install qt cmake
packaging/macos/fetch_indigo.sh
packaging/macos/build_imago.sh
cmake -S . -B build_mac -DCMAKE_PREFIX_PATH="$(brew --prefix qt)" -DCMAKE_BUILD_TYPE=Release
cmake --build build_mac -j
ctest --test-dir build_mac
packaging/macos/make_app.sh build_mac
packaging/macos/make_dmg.sh
```

- `fetch_indigo.sh` puts the Indigo 1.45.0 libraries, headers and licence in `indigo/` (gitignored). The libraries come from the `epam.indigo` arm64 wheel on PyPI.
- `build_imago.sh` builds Imago (image recognition, used by Insert Image) from source into
  `imago/` (gitignored). EPAM publishes no prebuilt Imago binaries for any platform, and Imago's
  own CMake build compiles its own vendored copies of Indigo and OpenCV from source too, so this
  is a slow step (expect it to take much longer than building QtSketch itself) and needs Xcode
  command line tools and git in addition to what's listed above. Imago's own upstream CI last
  built this exact pinned commit clean on a macOS runner, but that was not confirmed here against
  Apple silicon specifically -- build it once and confirm Insert Image works before relying on it.
- `make_app.sh` copies the build to `dist/QtSketch.app` and bundles Qt, Indigo and Imago, so the app runs without Homebrew.
- `make_dmg.sh` writes `dist/QtSketch-0.1-arm64.dmg` and its SHA-256.

## Limits

- arm64 only. The minimum macOS is whatever the bundled Qt was built for; Homebrew Qt on macOS 26 gives 26.0.
- The app is signed ad hoc, not notarised. Gatekeeper blocks it on first open; right-click the app and choose Open.
- The UI uses Windows fonts (Segoe UI, Segoe Fluent Icons, Consolas). On macOS Qt substitutes other fonts and some icons may be blank.
