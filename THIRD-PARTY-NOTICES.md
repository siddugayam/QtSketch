# Third-Party Notices

QtSketch uses the following third-party components.

## Qt (The Qt Company)

- **What is used:** Qt 6 Core, Gui, Qml, Quick, Quick Controls 2, Concurrent, Svg and
  SvgWidgets, linked dynamically.
- **Source:** https://www.qt.io and https://code.qt.io
- **License:** GNU Lesser General Public License v3.0 (open-source licensing of these modules)
- **Note:** Distributed builds ship the unmodified Qt libraries as separate files, so they can
  be replaced. The macOS app bundle also contains libraries that Homebrew's Qt depends on
  (for example ICU, OpenSSL, FreeType, HarfBuzz, GLib, libpng, zstd), each under its own
  licence.

## Indigo Toolkit (EPAM Systems)

- **What is used:** the Indigo C API (indigo, indigo-inchi, indigo-renderer) for layout,
  aromatisation, structure checks, CIP stereodescriptors, format conversion, rendering and
  biopolymer sequences. The binaries are not in the repository; they are placed in `indigo/`
  at build time.
- **Source:** https://github.com/epam/Indigo
- **License:** Apache License, Version 2.0
- **Copyright:** © EPAM Systems, Inc.

## InChI (IUPAC and InChI Trust)

- **What is used:** the InChI library (software v1.06) built into Indigo's `indigo-inchi`
  library, for InChI and InChIKey output.
- **Source:** the InChI source used by Indigo 1.45.0 is at
  https://github.com/epam/Indigo/tree/indigo-1.45.0/third_party/inchi
- **License:** IUPAC/InChI-Trust InChI Licence No. 1.0. QtSketch calls the library and does not
  modify it. Copies of QtSketch that include the `indigo-inchi` library must include this
  licence and point to the source above.
- **Copyright:** © IUPAC and InChI Trust Limited

## Imago (EPAM Systems)

- **What is used:** the Imago C API for recognising structures in images (Insert Image).
  Windows MinGW build only; the binaries are placed in `imago/` at build time.
- **Source:** https://github.com/epam/Imago
- **License:** Apache License, Version 2.0
- **Copyright:** © EPAM Systems, Inc.
- **Note:** Imago v1 was released under GPLv3 before the 2.1 relicense. QtSketch uses a build
  made from the current Apache-2.0 source, per the upstream `LICENSE` and `LICENSE-history`
  files. Imago builds on OpenCV, which is also Apache-2.0.

## Ketcher (EPAM Systems)

- **What is used:** the SVG icons in `icons/` come from the Ketcher user interface. The
  structure templates in `templates/` (`library.sdf`, `fg.sdf`, `salts-and-solvents.sdf`) are
  Ketcher's template files. `assets/monomer_library.ket` is the monomer library from the Indigo
  repository (`data/molecules/basic/`).
- **Source:** https://github.com/epam/ketcher
- **License:** Apache License, Version 2.0
- **Copyright:** © EPAM Systems, Inc.
- **Modifications:** `currentColor` fills and strokes were replaced with a fixed colour
  (#202020), because Qt's SVG renderer does not resolve CSS `currentColor`.

## Material Design Icons (Pictogrammers)

- **What is used:** the icon font `fonts/materialdesignicons-webfont.ttf`.
- **Source:** https://github.com/Templarian/MaterialDesign-Webfont
- **License:** Pictogrammers Free License; the font files are under the Apache License,
  Version 2.0.

## Libraries in the macOS app bundle

The macOS app (`packaging/macos/make_app.sh`) copies Homebrew's Qt and the libraries it depends
on into `QtSketch.app/Contents/Frameworks`. Each library is a separate, unmodified file. Licences
as declared by Homebrew:

| Library | License |
|---|---|
| Qt 6 frameworks and plugins | LGPL-3.0 (also GPL-2.0 / GPL-3.0 / commercial) |
| ICU | Unicode-3.0 |
| OpenSSL 3 | Apache-2.0 |
| FreeType | FreeType License (FTL) |
| HarfBuzz, Brotli, md4c | MIT |
| graphite2 | MIT, MPL-2.0, LGPL-2.1+ or GPL-2.0+ (choice) |
| GLib, GIO, GObject, libintl (gettext) | LGPL-2.1-or-later |
| D-Bus | AFL-2.1 or GPL-2.0-or-later (choice) |
| libpng | libpng-2.0 |
| libjpeg-turbo | IJG, Zlib and BSD-3-Clause |
| libtiff | libtiff |
| libwebp, libsharpyuv, PCRE2, double-conversion | BSD-3-Clause |
| zstd | BSD-3-Clause (or GPL-2.0) |
| liblzma (xz) | 0BSD |
| libb2 | CC0-1.0 |

The bundle carries the licence texts for Apache-2.0, LGPL-3.0, GPL-3.0, LGPL-2.1 and the InChI
Licence in
`Contents/Resources/licenses/`.

## Apache License 2.0

Indigo, Imago, Ketcher and the Material Design Icons font are licensed under the Apache
License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0. Copies of the licence text
come with the Indigo and Imago distributions (`indigo/LICENSE`, `imago/LICENSE`) placed at
build time.

The names of these projects and of EPAM Systems are not used to endorse or promote QtSketch.
