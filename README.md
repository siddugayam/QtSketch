# ![QtSketch](icons/app_icon_48.png) QtSketch

**Fast, native desktop editor for chemical structures, reactions, and biopolymers**

[![Download macOS DMG](https://img.shields.io/badge/Download%20macOS-Apple%20Silicon%20(.dmg%20%7C%2051%20MB)-000000?logo=apple&logoColor=white)](releases/QtSketch-0.1-arm64.dmg)
[![Linux AppImage](https://img.shields.io/badge/Linux-x86__64%20%7C%20arm64%20(.AppImage)-FCC624?logo=linux&logoColor=black)](#-quick-setup-and-build-macos-and-linux)
[![Windows x64](https://img.shields.io/badge/Windows-x64%20(MinGW%20%7C%20MSVC)-0078D6?logo=windows&logoColor=white)](#-building-on-windows-mingw)
[![Qt 6](https://img.shields.io/badge/Qt-6.6%2B%20%7C%20QML-41CD52?logo=qt&logoColor=white)](https://www.qt.io/)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-00599C?logo=cplusplus&logoColor=white)](https://isocpp.org/)
[![Indigo](https://img.shields.io/badge/Indigo-1.45.0-007ACC)](https://github.com/epam/Indigo)
[![Imago OCR](https://img.shields.io/badge/Imago-Structure%20OCR-5C2D91)](https://github.com/epam/Imago)
[![CTest](https://img.shields.io/badge/CTest-28%2F28%20passing-brightgreen)](#-tests)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

---

**QtSketch** is a cross-platform desktop chemical structure sketcher built with **Qt 6 (QML)** and **C++20** on the [EPAM Indigo](https://github.com/epam/Indigo) cheminformatics toolkit. Draw and edit small molecules, multi-step reaction schemes, and biopolymer sequences, recognize chemical structures directly from images, and generate systematic IUPAC names offline.

## 📥 Direct Downloads

| Installer | Architecture | Format & Size | SHA-256 Checksum |
| :--- | :--- | :--- | :--- |
| 🍎 **[QtSketch-0.1-arm64.dmg](releases/QtSketch-0.1-arm64.dmg)** | macOS Apple Silicon (`arm64`) | Self-contained `.dmg` (`51 MB`, LZMA) | [`QtSketch-0.1-arm64.dmg.sha256`](releases/QtSketch-0.1-arm64.dmg.sha256) |
| 🐧 **Linux `.AppImage`** | `x86_64` & `arm64` | Built via `./install_and_build.sh` or [GitHub Releases](../../releases) | Generated `.sha256` |

---

## ✨ Key Features

| Feature | Capabilities |
| :--- | :--- |
| 🧪 **2D Structure Drawing** | Atoms, single/double/triple/wedge/hash/wavy/dative/aromatic bonds, rings, ring templates, functional-group abbreviations, brackets, and Markush R-groups |
| 🧬 **Reactions & Biopolymers** | Multi-tail & curved reaction arrows, atom-atom mapping, and full biopolymer sequence editing (peptides, RNA, DNA, HELM, and FASTA) |
| 🔬 **Cheminformatics & Stereo** | 2D clean-up/layout, aromatisation, valence & structure checks, CIP *(R/S, E/Z)* stereodescriptors, SMILES, SMARTS search, InChI, Molfile, KET, and CML |
| 🔍 **Image OCR (Imago)** | Optical chemical structure recognition (*Insert Image*) powered by [EPAM Imago](https://github.com/epam/Imago) |
| 🏷️ **IUPAC Name Generator** | Built-in **IUPAC Blue Book** systematic name generator (`src/app/IupacNamer.cpp`) verified against a golden test set |
| 🖼️ **Publication Export** | Crisp vector and high-DPI raster export to **SVG**, **PDF**, and **PNG** (with ACS/RSC style sheets), plus batch export to **SDF**, **RDF**, **SMILES**, and **CML** |

---

## 🖥️ Supported Platforms

| Platform | Toolchain | Indigo | Imago OCR | Status |
| :--- | :--- | :--- | :--- | :--- |
| 🍎 **macOS `arm64`** | Homebrew Qt 6, Xcode CLT | EPAM `epam.indigo` 1.45.0 wheel | ✅ `packaging/macos/build_imago.sh` | Prebuilt `.dmg` in [`releases/`](releases/), native macOS window chrome, `28/28` tests pass |
| 🐧 **Linux `x86_64`, `arm64`** | System Qt 6.6+ or Qt 6.8 (`aqtinstall`), native or Docker | EPAM `epam.indigo` 1.45.0 wheel | ✅ `packaging/linux/build_imago.sh` | One-step `./install_and_build.sh`, Docker `.AppImage`, and GitHub Actions CI |
| 🪟 **Windows `x64`** *(MinGW)* | Qt 6.6+ MinGW 64-bit kit | MinGW build in `indigo/` | ✅ `imago/lib/imago.dll` | Primary Windows development & release toolchain |
| 🪟 **Windows `x64`** *(MSVC)* | Qt 6.6+ MSVC kit | EPAM `indigo-libs-windows-x86_64` | ➖ Optional (`QTSKETCH_WITH_IMAGO=OFF`) | CMake MSVC configuration supported |

---

## ⚡ Quick Setup and Build (macOS and Linux)

[`install_and_build.sh`](install_and_build.sh) installs missing OS packages, fetches Indigo (with SHA-256 verification), compiles Imago from source, builds QtSketch, runs all 28 CTest suites, and packages a self-contained `.dmg` (macOS) or `.AppImage` (Linux):

```bash
./install_and_build.sh
```

Useful flags:

```bash
./install_and_build.sh --skip-imago     # Fast build without Imago OCR (~30 s)
./install_and_build.sh --skip-package   # Compile and test only; skip .dmg / .AppImage
./install_and_build.sh --no-deps        # Skip brew / apt / dnf / pacman / zypper install
./install_and_build.sh --run            # Launch QtSketch immediately after building
```

> [!TIP]
> On Linux distributions whose system Qt or CMake is older than required (such as Ubuntu 22.04 LTS and Ubuntu 24.04 LTS), `install_and_build.sh` automatically bootstraps CMake ≥ 3.24 and Qt 6.8.3 locally inside `.qt/` via `aqtinstall` without modifying system packages.

Prebuilt macOS `.dmg` and Linux `x86_64` / `aarch64` `.AppImage` installers are also built automatically by [`.github/workflows/release.yml`](.github/workflows/release.yml) on `workflow_dispatch` and published to **GitHub Releases** whenever a `v*` tag is pushed.

---

## 📦 Requirements (Manual Builds)

- **Qt 6.6+** with `Core`, `Gui`, `Qml`, `Quick`, `QuickControls2`, `Concurrent`, `Svg`, `SvgWidgets`, and `Test`
- **CMake 3.24+** and **Ninja** (or MinGW Makefiles on Windows)
- **C++20 compiler** matching your Qt kit (Apple Clang on macOS, GCC/Clang on Linux, MinGW GCC or MSVC on Windows)
- **Indigo 1.43+** C API headers and shared libraries in `indigo/` (fetched automatically on macOS/Linux)
- **Imago** *(optional, enabled automatically when `imago/lib/` is present)* — built from source on macOS/Linux via `packaging/{macos,linux}/build_imago.sh`

---

## 🛠️ Platform Build Guides

### 🍎 Building on macOS (Manual)

```bash
brew install qt cmake ninja
packaging/macos/fetch_indigo.sh
packaging/macos/build_imago.sh          # optional; omit for a faster build without OCR
cmake -S . -B build_mac -G Ninja -DCMAKE_PREFIX_PATH="$(brew --prefix qt)" -DCMAKE_BUILD_TYPE=Release
cmake --build build_mac
QT_QPA_PLATFORM=offscreen ctest --test-dir build_mac --output-on-failure
packaging/macos/make_app.sh build_mac
packaging/macos/make_dmg.sh
```

Produces `dist/QtSketch.app` and `dist/QtSketch-0.1-arm64.dmg` (plus `.sha256`).

> [!NOTE]
> `QtSketch.app` is ad hoc code-signed rather than Apple-notarised. On first open after mounting the `.dmg`, right-click **QtSketch.app → Open** or run `xattr -cr /Applications/QtSketch.app`. See [packaging/macos/README.md](packaging/macos/README.md) for details.

### 🐧 Building on Linux (Manual or Docker)

Build a portable `x86_64` or `arm64` `.AppImage` inside Docker from any host:

```bash
packaging/linux/build_in_docker.sh amd64     # or arm64
```

Or build natively on Ubuntu/Debian, Fedora/RHEL, Arch, or openSUSE:

```bash
./install_and_build.sh
```

See [packaging/linux/README.md](packaging/linux/README.md) for distribution package lists and AppImage details.

### 🪟 Building on Windows (MinGW)

#### 1. Install Windows tools

- **Qt 6.6+** with the **MinGW 64-bit** kit from the [Qt online installer](https://www.qt.io/download-qt-installer) or `aqtinstall`.
- Use the **exact MinGW toolchain** bundled with that Qt kit (for example `C:\Qt\Tools\mingw1310_64`). Do not mix MinGW from Git Bash, Conda, or MSYS2 — mismatched `libstdc++` DLLs cause `STATUS_ENTRYPOINT_NOT_FOUND` at start-up.
- **CMake 3.24+** (`winget install Kitware.CMake`).

#### 2. Add Indigo and Imago

Place MinGW-compatible builds of Indigo and Imago in the repository root:

```text
indigo/
    LICENSE
    api/c/indigo/            indigo.h
    api/c/indigo-inchi/
    api/c/indigo-renderer/
    lib/indigo.dll, indigo-inchi.dll, indigo-renderer.dll
imago/
    LICENSE
    include/imago_c.h        (copy from packaging/imago_c.h)
    lib/imago.dll
```

#### 3. Configure and build on Windows

From the repository root in PowerShell:

```powershell
cmake -B build -G "MinGW Makefiles" `
    -DCMAKE_PREFIX_PATH="C:\Qt\6.11.1\mingw_64" `
    -DCMAKE_C_COMPILER="C:\mingw64\bin\gcc.exe" `
    -DCMAKE_CXX_COMPILER="C:\mingw64\bin\g++.exe"
$env:PATH = "C:\mingw64\bin;$env:PATH"
mingw32-make -C build -j4
```

#### 4. Launch on Windows

```bat
@echo off
set PATH=C:\Qt\6.11.1\mingw_64\bin;%~dp0indigo\lib;%PATH%
start "" "%~dp0build\sketch.exe"
```

If you switch generators or toolchains, delete `build\CMakeCache.txt` and `build\CMakeFiles` before reconfiguring.

---

## 🧪 Tests

**28 automated CTest suites** are registered with CMake and built alongside the app:

```bash
QT_QPA_PLATFORM=offscreen ctest --test-dir build --output-on-failure
```

| Area | Test Suites |
| :--- | :--- |
| 🏷️ **IUPAC Naming** | `iupac_namer_test`, `iupac_golden_test`, `naming_validation_test`, `fused_ring_orientation_test`, `fused_ring_direction_detector_test`, `fused_system_topology_test`, `peripheral_numbering_test`, `system_locant_test`, `fusion_construction_test`, `retained_ring_registry_test`, `retained_ring_matcher_test` |
| 🧬 **Document & Biopolymers** | `editable_molecule_test`, `document_state_test`, `biopolymer_state_test`, `placement_engines_test`, `bond_angle_suggester_test`, `clipboard_preview_test` |
| 💾 **Files & Settings** | `project_file_test`, `settings_service_test`, `sdf_batch_test` |
| 🎨 **Rendering & Export** | `render_scene_test`, `export_manager_test`, `render_primitives_test`, `render_primitives_to_variant_test`, `biopolymer_sequence_view_test` |
| 🔎 **Guards & Image OCR** | `bio_input_guards_test`, `imago_service_test`, `imago_service_stub_test` |

---

## 🗂️ Project Layout

```text
QtSketch/
├── src/app/               C++20 backend: DocumentState, EditableMolecule, IndigoService,
│                          ImagoService, IupacNamer, vector/raster rendering & export
├── src/v8_process.*       QML command bridge (historical name; invokes DocumentState directly)
├── *.qml, components/,    Qt Quick / QML user interface, toolbars, canvas & dialogs
│   dialogs/, js/
├── templates/             Chemical ring & functional-group template library
├── assets/, icons/,       Monomer library (monomer_library.ket), 250+ SVG icons, icon font
│   fonts/
├── tests/                 28 CTest executables and IUPAC Blue Book golden datasets
├── install_and_build.sh   One-command OS setup, build, test, and .dmg/.AppImage packager
├── releases/              Prebuilt macOS Apple Silicon .dmg installer and SHA-256 checksum
├── packaging/
│   ├── macos/             Indigo fetcher, Imago source builder, .app & .dmg scripts
│   └── linux/             Indigo fetcher, Imago source builder, Docker & AppImage scripts
├── .github/workflows/     Automated macOS (.dmg) & Linux (.AppImage) CI and GitHub Releases
└── docs/iupac/            Architecture and design notes for the IUPAC name generator
```

---

## ⚖️ Licence

QtSketch is licensed under the **[Apache License, Version 2.0](LICENSE)** — see [NOTICE](NOTICE).
Third-party components and licences are listed in [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).
