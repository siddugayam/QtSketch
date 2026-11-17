#!/bin/bash
# Makes a self-contained, ad hoc signed dist/QtSketch.app from a CMake build of the repo.
# Usage: packaging/macos/make_app.sh [build-dir]   (default: build)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD="${1:-$ROOT/build}"
APP="$ROOT/dist/QtSketch.app"
MACDEPLOYQT="${MACDEPLOYQT:-$(command -v macdeployqt || echo /opt/homebrew/opt/qt/bin/macdeployqt)}"

[ -d "$BUILD/QtSketch.app" ] || { echo "no $BUILD/QtSketch.app; build first"; exit 1; }
[ -x "$MACDEPLOYQT" ] || { echo "macdeployqt not found; set MACDEPLOYQT"; exit 1; }

rm -rf "$APP"; mkdir -p "$ROOT/dist"
cp -R "$BUILD/QtSketch.app" "$APP"
BIN="$APP/Contents/MacOS/QtSketch"

mkdir -p "$APP/Contents/Frameworks"
cp "$ROOT"/indigo/lib/libindigo*.dylib "$APP/Contents/Frameworks/"

# IndigoService looks for <exe dir>/../indigo and DocumentManager walks up to templates/.
# Both are found under Contents/, which is inside the code signature.
mkdir -p "$APP/Contents/indigo/data/molecules/basic"
cp "$ROOT/assets/monomer_library.ket" "$APP/Contents/indigo/data/molecules/basic/"
cp -R "$ROOT/templates" "$APP/Contents/templates"
mkdir -p "$APP/Contents/Resources/licenses"
cp "$ROOT/NOTICE" "$ROOT/THIRD-PARTY-NOTICES.md" "$APP/Contents/Resources/licenses/"
cp "$ROOT/indigo/LICENSE" "$APP/Contents/Resources/licenses/Apache-2.0.txt"
cp "$ROOT/indigo/LICENCE-InChI" "$APP/Contents/Resources/licenses/InChI-Licence-1.0.txt"
# Homebrew ships no licence texts. Qt is LGPL-3.0 (which builds on GPL-3.0); GLib and libintl are LGPL-2.1.
for l in LGPL-3.0-only GPL-3.0-only LGPL-2.1-only; do
    curl -fsSL --max-time 60 "https://raw.githubusercontent.com/spdx/license-list-data/main/text/$l.txt" \
        -o "$APP/Contents/Resources/licenses/$l.txt"
done

# macdeployqt can exit non-zero on harmless warnings, so check its result instead.
"$MACDEPLOYQT" "$APP" -qmldir="$ROOT" -no-strip >/dev/null 2>&1 || true
for fw in QtCore QtQuick; do
    [ -d "$APP/Contents/Frameworks/$fw.framework" ] || { echo "macdeployqt did not bundle $fw"; exit 1; }
done

# Homebrew Qt: macdeployqt leaves some copied libraries pointing at /opt/homebrew.
find "$APP/Contents" -type f \( -name "*.dylib" -o -path "*.framework/Versions/A/*" -o -path "*MacOS/*" \) | while read -r f; do
    file "$f" | grep -q "Mach-O" || continue
    chmod u+w "$f"
    otool -L "$f" | awk 'NR>1{print $1}' | { grep "^/opt/homebrew.*\.framework/" || true; } | while read -r dep; do
        install_name_tool -change "$dep" "@executable_path/../Frameworks/${dep##*/lib/}" "$f" 2>/dev/null || true
    done
done

# Drop absolute build rpaths so only the bundled Qt and Indigo are loaded.
otool -l "$BIN" | awk '/LC_RPATH/{getline; getline; print $2}' | { grep '^/' || true; } | while read -r rp; do
    install_name_tool -delete_rpath "$rp" "$BIN"
done
install_name_tool -add_rpath @executable_path/../Frameworks "$BIN" 2>/dev/null || true
rm -f "$APP/Contents/PlugIns/imageformats/libqpdf.dylib"
printf '[Paths]\nPlugins = PlugIns\n' > "$APP/Contents/Resources/qt.conf"

# Minimum macOS = highest minos among the bundled binaries.
MINOS=$(find "$APP/Contents" -type f \( -name "*.dylib" -o -path "*.framework/Versions/A/*" -o -path "*MacOS/*" \) \
    -exec otool -l {} \; 2>/dev/null | awk '/minos/{print $2}' | sort -V | tail -1)
PL="$APP/Contents/Info.plist"
chmod u+w "$PL"
/usr/libexec/PlistBuddy -c "Set :LSMinimumSystemVersion $MINOS" "$PL" 2>/dev/null \
    || /usr/libexec/PlistBuddy -c "Add :LSMinimumSystemVersion string $MINOS" "$PL"
/usr/libexec/PlistBuddy -c "Add :NSHighResolutionCapable bool true" "$PL" 2>/dev/null || true

ICS="$(mktemp -d)/QtSketch.iconset"; mkdir -p "$ICS"
for sz in 16 32 128 256 512; do
    sips -z $sz $sz "$ROOT/icons/app_icon.png" --out "$ICS/icon_${sz}x${sz}.png" >/dev/null
    sips -z $((sz*2)) $((sz*2)) "$ROOT/icons/app_icon.png" --out "$ICS/icon_${sz}x${sz}@2x.png" >/dev/null
done
iconutil -c icns "$ICS" -o "$APP/Contents/Resources/QtSketch.icns"
/usr/libexec/PlistBuddy -c "Set :CFBundleIconFile QtSketch" "$PL" 2>/dev/null \
    || /usr/libexec/PlistBuddy -c "Add :CFBundleIconFile string QtSketch" "$PL"

codesign --force --deep --sign - "$APP"
codesign --verify --deep --strict "$APP"
echo "built $APP (macOS $MINOS or later)"
