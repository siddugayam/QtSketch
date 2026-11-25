#!/bin/bash
# Builds a self-contained, ad hoc signed dist/QtSketch.app from a CMake build of the repo.
# Usage: packaging/macos/make_app.sh [build-dir]   (default: build_mac, fallback: build)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
if [ $# -ge 1 ]; then
    BUILD="$1"
elif [ -d "$ROOT/build_mac/QtSketch.app" ]; then
    BUILD="$ROOT/build_mac"
else
    BUILD="$ROOT/build"
fi
case "$BUILD" in
    /*) ;;
    *) BUILD="$ROOT/$BUILD" ;;
esac

APP="$ROOT/dist/QtSketch.app"
QT_PREFIX="${QT_PREFIX:-$(brew --prefix qt 2>/dev/null || echo /opt/homebrew/opt/qt)}"
MACDEPLOYQT="${MACDEPLOYQT:-$(command -v macdeployqt || echo "$QT_PREFIX/bin/macdeployqt")}"

[ -d "$BUILD/QtSketch.app" ] || { echo "No $BUILD/QtSketch.app found; build first."; exit 1; }
[ -x "$MACDEPLOYQT" ] || { echo "macdeployqt not found; set MACDEPLOYQT or install Homebrew qt."; exit 1; }

rm -rf "$APP"
mkdir -p "$ROOT/dist"
cp -R "$BUILD/QtSketch.app" "$APP"
BIN="$APP/Contents/MacOS/QtSketch"

mkdir -p "$APP/Contents/Frameworks"
cp "$ROOT"/indigo/lib/libindigo*.dylib "$APP/Contents/Frameworks/"
if compgen -G "$ROOT/imago/lib/libimago*.dylib" >/dev/null; then
    cp "$ROOT"/imago/lib/libimago*.dylib "$APP/Contents/Frameworks/"
fi

for dylib in "$APP/Contents/Frameworks"/libindigo*.dylib "$APP/Contents/Frameworks"/libimago*.dylib; do
    [ -f "$dylib" ] || continue
    chmod u+w "$dylib"
    install_name_tool -id "@rpath/$(basename "$dylib")" "$dylib" 2>/dev/null || true
done

# Bundle runtime data inside Contents/Resources (required by Apple's CFBundle code-signing seal).
mkdir -p "$APP/Contents/Resources/indigo/data/molecules/basic"
cp "$ROOT/assets/monomer_library.ket" "$APP/Contents/Resources/indigo/data/molecules/basic/"
cp -R "$ROOT/templates" "$APP/Contents/Resources/templates"

LIC="$APP/Contents/Resources/licenses"
mkdir -p "$LIC"
cp "$ROOT/LICENSE" "$ROOT/NOTICE" "$ROOT/THIRD-PARTY-NOTICES.md" "$LIC/"
cp "$ROOT/indigo/LICENSE" "$LIC/Indigo-Apache-2.0.txt"
cp "$ROOT/indigo/LICENCE-InChI" "$LIC/InChI-Licence-1.0.txt"
if compgen -G "$ROOT/imago/lib/libimago*.dylib" >/dev/null; then
    cp "$ROOT/imago/LICENSE" "$LIC/Imago-Apache-2.0.txt"
    cp "$ROOT/imago/LICENSE-OpenCV" "$LIC/OpenCV-Apache-2.0.txt"
fi
SPDX_CACHE="$ROOT/dist/.spdx-cache"
mkdir -p "$SPDX_CACHE"
for l in LGPL-3.0-only GPL-3.0-only LGPL-2.1-only; do
    if [ ! -s "$SPDX_CACHE/$l.txt" ]; then
        curl -fsSL --retry 3 --max-time 60 \
            "https://raw.githubusercontent.com/spdx/license-list-data/main/text/$l.txt" \
            -o "$SPDX_CACHE/$l.txt"
    fi
    cp "$SPDX_CACHE/$l.txt" "$LIC/$l.txt"
    [ -s "$LIC/$l.txt" ] || { echo "Missing required license file: $l.txt"; exit 1; }
done

# macdeployqt can return non-zero on harmless warnings; verify bundled frameworks explicitly.
"$MACDEPLOYQT" "$APP" -qmldir="$ROOT" -no-strip >/dev/null 2>&1 || true
for fw in QtCore QtGui QtQml QtQuick; do
    [ -d "$APP/Contents/Frameworks/$fw.framework" ] || { echo "macdeployqt did not bundle $fw"; exit 1; }
done

# Fix any remaining Homebrew or temporary build paths in Mach-O load commands (multi-pass for newly copied dylibs).
for _pass in 1 2 3; do
    while IFS= read -r f; do
        file "$f" | grep -q "Mach-O" || continue
        chmod u+w "$f"
        case "$f" in
            "$APP"/Contents/Frameworks/*.framework/Versions/*/*)
                fw_rel="${f#"$APP/Contents/Frameworks/"}"
                install_name_tool -id "@rpath/$fw_rel" "$f" 2>/dev/null || true
                ;;
            "$APP"/Contents/Frameworks/*.dylib)
                install_name_tool -id "@rpath/$(basename "$f")" "$f" 2>/dev/null || true
                ;;
        esac
        while IFS= read -r dep; do
            [ -n "$dep" ] || continue
            case "$dep" in
                /opt/homebrew/*.framework/*|/usr/local/*.framework/*)
                    rel="${dep##*/lib/}"
                    install_name_tool -change "$dep" "@executable_path/../Frameworks/$rel" "$f" 2>/dev/null || true
                    ;;
                /opt/homebrew/*.dylib|/usr/local/*.dylib)
                    base="$(basename "$dep")"
                    if [ ! -f "$APP/Contents/Frameworks/$base" ] && [ -f "$dep" ]; then
                        cp -L "$dep" "$APP/Contents/Frameworks/$base"
                        chmod u+w "$APP/Contents/Frameworks/$base"
                        install_name_tool -id "@rpath/$base" "$APP/Contents/Frameworks/$base" 2>/dev/null || true
                    fi
                    install_name_tool -change "$dep" "@executable_path/../Frameworks/$base" "$f" 2>/dev/null || true
                    ;;
                */libimago*.dylib|*/libindigo*.dylib)
                    case "$dep" in
                        @*) ;;
                        *) install_name_tool -change "$dep" "@rpath/$(basename "$dep")" "$f" 2>/dev/null || true ;;
                    esac
                    ;;
            esac
        done < <(otool -L "$f" | awk 'NR>1{print $1}')
    done < <(find "$APP/Contents" -type f \( -name "*.dylib" -o -path "*.framework/Versions/*/*" -o -path "*/MacOS/*" \))
done

# Drop absolute build RPATHs so only bundled Frameworks are searched at runtime.
otool -l "$BIN" | awk '/LC_RPATH/{getline; getline; print $2}' | { grep '^/' || true; } | while read -r rp; do
    [ -n "$rp" ] && install_name_tool -delete_rpath "$rp" "$BIN" 2>/dev/null || true
done
install_name_tool -add_rpath @executable_path/../Frameworks "$BIN" 2>/dev/null || true

rm -f "$APP/Contents/PlugIns/imageformats/libqpdf.dylib"
# Prefix = .. sets the bundle root to QtSketch.app/Contents so PlugIns and
# Resources/qml resolve properly when launched standalone without Homebrew Qt.
printf '[Paths]\nPrefix = ..\nPlugins = PlugIns\nImports = Resources/qml\nQmlImports = Resources/qml\nQml2Imports = Resources/qml\n' > "$APP/Contents/Resources/qt.conf"

MINOS=$(find "$APP/Contents" -type f \( -name "*.dylib" -o -path "*.framework/Versions/A/*" -o -path "*/MacOS/*" \) \
    -exec otool -l {} \; 2>/dev/null | awk '/minos/{print $2}' | sort -V | tail -1)
MINOS="${MINOS:-14.0}"
PL="$APP/Contents/Info.plist"
chmod u+w "$PL"
/usr/libexec/PlistBuddy -c "Set :LSMinimumSystemVersion $MINOS" "$PL" 2>/dev/null \
    || /usr/libexec/PlistBuddy -c "Add :LSMinimumSystemVersion string $MINOS" "$PL"
/usr/libexec/PlistBuddy -c "Add :NSHighResolutionCapable bool true" "$PL" 2>/dev/null || true

TMP_ICON="$(mktemp -d)"
ICS="$TMP_ICON/QtSketch.iconset"
mkdir -p "$ICS"
for sz in 16 32 128 256 512; do
    sips -z "$sz" "$sz" "$ROOT/icons/app_icon.png" --out "$ICS/icon_${sz}x${sz}.png" >/dev/null
    sips -z "$((sz*2))" "$((sz*2))" "$ROOT/icons/app_icon.png" --out "$ICS/icon_${sz}x${sz}@2x.png" >/dev/null
done
iconutil -c icns "$ICS" -o "$APP/Contents/Resources/QtSketch.icns" 2>/dev/null \
    || sips -s format icns -z 512 512 "$ROOT/icons/app_icon.png" --out "$APP/Contents/Resources/QtSketch.icns" >/dev/null
rm -rf "$TMP_ICON"
/usr/libexec/PlistBuddy -c "Set :CFBundleIconFile QtSketch" "$PL" 2>/dev/null \
    || /usr/libexec/PlistBuddy -c "Add :CFBundleIconFile string QtSketch" "$PL"

codesign --force --deep --sign - "$APP"
codesign --verify --deep --strict "$APP"
echo "Built $APP (macOS $MINOS or later)"
