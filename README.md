# QtSketch

QtSketch is a desktop editor for chemical structures, built with Qt 6 (QML) and C++20 on the
[Indigo](https://github.com/epam/Indigo) cheminformatics toolkit. You can use it to draw and edit
molecules, reactions and biopolymer sequences. It also includes an IUPAC name generator written
from scratch.

## Features

- 2D drawing and editing of atoms, bonds, rings, templates, abbreviation groups and R-groups
- Reaction schemes and biopolymer sequences (peptide, RNA, DNA, HELM, FASTA)
- Layout and clean-up, aromatisation, structure checks, CIP stereodescriptors, SMILES, InChI,
  molfile, KET and CML through Indigo
- Export to PNG, SVG and PDF, plus batch export to SDF, RDF, SMILES and CML
- Recognition of structures in images (Insert Image), through
  [Imago](https://github.com/epam/Imago) (Windows MinGW build only)
- IUPAC name generation (`src/app/IupacNamer.cpp`), written against the IUPAC Blue Book and
  checked by a golden test set

## Supported platforms

| Platform | Toolchain | Indigo | Imago | Status |
|---|---|---|---|---|
| Windows x64 | Qt 6.6+ MinGW kit | MinGW build in `indigo/` | yes | Main development platform |
| Windows x64 | Qt 6.6+ MSVC kit | EPAM `indigo-libs-windows-x86_64` | no | CMake branch present |
| macOS arm64 | Homebrew Qt 6.12, Xcode CLT | EPAM `epam.indigo` 1.45.0 wheel | no | Builds, 21/21 tests pass |
| Linux | — | — | — | Not supported yet |

## Requirements

- Qt 6.6 or newer, with Core, Gui, Qml, Quick, QuickControls2, Concurrent, Svg, SvgWidgets and
  Test
- CMake 3.24 or newer
- A C++20 compiler matching the Qt kit (MinGW GCC or MSVC on Windows, Apple Clang on macOS)
- Indigo 1.43 or newer (C API headers and libraries), placed in the repository as described
  below. Indigo is not committed to git.
- Imago (optional, Windows MinGW only), used by Insert Image

## Building on Windows (MinGW)

### 1. Install the tools

- Qt 6.6 or newer with the MinGW 64-bit kit, from the
  [Qt online installer](https://www.qt.io/download-qt-installer) or `aqtinstall`.
- The MinGW toolchain that ships with that kit (for example `C:\Qt\Tools\mingw1310_64`).
  Do not use a different MinGW from Git Bash, Conda or MSYS2. A mismatched `libstdc++` causes
  `STATUS_ENTRYPOINT_NOT_FOUND` crashes at start-up.
- CMake 3.24 or newer (`winget install Kitware.CMake`).

### 2. Add Indigo and Imago

Both are Apache-2.0 projects by EPAM Systems. Build them from source with the same MinGW
toolchain, or copy a MinGW-compatible build, into this layout at the repository root:

```
indigo/
    LICENSE
    api/c/indigo/            indigo.h
    api/c/indigo-inchi/
    api/c/indigo-renderer/
    lib/indigo.dll, indigo-inchi.dll, indigo-renderer.dll
imago/
    LICENSE
    include/imago_c.h
    lib/imago.dll
```

### 3. Configure and build

From the repository root in PowerShell:

```powershell
cmake -B build -G "MinGW Makefiles" `
    -DCMAKE_PREFIX_PATH="C:\Qt\6.11.1\mingw_64" `
    -DCMAKE_C_COMPILER="C:\mingw64\bin\gcc.exe" `
    -DCMAKE_CXX_COMPILER="C:\mingw64\bin\g++.exe"
$env:PATH = "C:\mingw64\bin;$env:PATH"
mingw32-make -C build -j4
```

Change the paths to match your Qt kit and MinGW install. Give the generator explicitly.
Without it, CMake may pick a Visual Studio generator and then fail to find Qt.

The build copies the Indigo and Imago DLLs next to `sketch.exe`, and copies
`assets/monomer_library.ket` to `indigo/data/molecules/basic/`.

### 4. Run

Put Qt's `bin` folder and `indigo\lib` on `PATH`, then start `build\sketch.exe`. A small
launcher script (not committed) does this:

```bat
@echo off
set PATH=C:\Qt\6.11.1\mingw_64\bin;%~dp0indigo\lib;%PATH%
start "" "%~dp0build\sketch.exe"
```

If you change generator or toolchain, delete `build\CMakeCache.txt` and `build\CMakeFiles`
before configuring again.

## Building on macOS

```bash
brew install qt cmake
packaging/macos/fetch_indigo.sh
cmake -S . -B build_mac -DCMAKE_PREFIX_PATH="$(brew --prefix qt)" -DCMAKE_BUILD_TYPE=Release
cmake --build build_mac -j
ctest --test-dir build_mac
packaging/macos/make_app.sh build_mac
packaging/macos/make_dmg.sh
```

This produces `dist/QtSketch.app` and `dist/QtSketch-0.1-arm64.dmg`. See
[packaging/macos/README.md](packaging/macos/README.md) for details and current limits
(ad hoc signature, no Imago, Windows fonts).

## Tests

21 Qt Test executables are registered with CTest and built with the app:

```
ctest --test-dir build --output-on-failure
```

| Area | Tests |
|---|---|
| IUPAC naming | `iupac_namer_test`, `iupac_golden_test`, `naming_validation_test`, `fused_ring_orientation_test`, `fused_ring_direction_detector_test`, `fused_system_topology_test`, `retained_ring_registry_test`, `retained_ring_matcher_test` |
| Document model | `editable_molecule_test`, `document_state_test`, `placement_engines_test`, `bond_angle_suggester_test`, `clipboard_preview_test` |
| Files and settings | `project_file_test`, `settings_service_test`, `sdf_batch_test` |
| Rendering and export | `render_scene_test`, `export_manager_test`, `render_primitives_test`, `render_primitives_to_variant_test`, `biopolymer_sequence_view_test` |

## Project layout

```
src/app/            C++ backend: DocumentState, EditableMolecule, IndigoService,
                    ImagoService, IupacNamer, rendering and export
src/v8_process.*    Command surface exposed to QML (the name is historical; it calls
                    DocumentState directly)
*.qml, components/, dialogs/   User interface
js/                 Small QML helper scripts
tests/              Qt Test suites and IUPAC golden data
templates/          Structure template library
assets/, icons/, fonts/        Monomer library, toolbar icons, icon font
packaging/macos/    macOS build and packaging scripts
docs/iupac/         Notes on the IUPAC namer
```

An earlier version ran its chemistry through an embedded QuickJS script engine. That engine is
retired and is no longer in the repository.

## Licence

The repository does not yet state a licence for QtSketch itself. Third-party components and
their licences are listed in [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).
