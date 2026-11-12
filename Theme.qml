pragma Singleton
import QtQuick

QtObject {
    // Toggled from the status bar; persisted via MainWindow's ui Settings.
    property bool compactDensity: true
    property bool cpkColorsEnabled: true

    // ── Fluent (Windows 11) light palette ──────────────────────────────────
    // Names on the right are the WinUI 3 resource keys these mirror, so any
    // value here can be checked against the Fluent design token it comes from.
    readonly property color surface: "#FFFFFF"            // LayerFillColorDefault
    readonly property color background: "#F3F3F3"         // SolidBackgroundFillColorBase
    readonly property color accent: "#0067C0"             // AccentFillColorDefault
    readonly property color accentLight: "#EFF6FC"        // AccentFillColorSelectedBackground
    readonly property color outline: "#E5E5E5"            // CardStrokeColorDefault
    readonly property color textPrimary: "#1B1B1B"        // TextFillColorPrimary
    readonly property color textSecondary: "#5D5D5D"      // TextFillColorSecondary
    readonly property color textTertiary: "#767676"       // TextFillColorTertiary
    readonly property color rgroupColor: "#8764B8"        // Fluent purple
    readonly property color hover: "#F5F5F5"              // ControlFillColorSecondary
    readonly property color pressed: "#EDEDED"            // ControlFillColorTertiary
    readonly property color selected: "#EFF6FC"           // SubtleFillColorSelected
    // Ruler numerals are 12px text on `background`; TextFillColorSecondary
    // clears 4.5:1 there, the old warm gray (#9E9E8E) measured 2.40:1.
    readonly property color rulerColor: "#5D5D5D"
    readonly property color workspaceBackground: "#EAEAEA" // canvas surround, one step under base
    readonly property color error: "#C42B1C"              // SystemFillColorCritical
    readonly property color success: "#0F7B0F"            // SystemFillColorSuccess
    readonly property color caution: "#9D5D00"            // SystemFillColorCaution
    readonly property color closeButtonHover: "#C42B1C"   // Fluent caption close
    readonly property color smokeFill: "#4D000000"        // SmokeFillColorDefault (30% black)
    readonly property color badgeAam: "#B7540A"
    readonly property color badgeReactingCenter: "#8764B8"
    readonly property color badgeText: "#FFFFFF"

    // ── Typography — Windows 11 type ramp ──────────────────────────────────
    // Segoe UI Variable is the Windows 11 UI face; Text and Display are its
    // optical sizes (Text is hinted for small sizes, Display for headings).
    // Both fall back to Segoe UI on Windows 10, which is metrically close.
    // Cascadia Mono is the Fluent monospace face, used for property values,
    // SMILES strings and the status readout where digits must align.
    // One family name per property, never a comma list: QML's font.family takes
    // a single family and treats "A, B" as a literal name that matches nothing,
    // silently falling back. Qt's own font matching handles absent families.
    readonly property string fontFamily: "Segoe UI Variable Text"
    readonly property string fontDisplay: "Segoe UI Variable Display"
    // Consolas rather than Fluent's Cascadia Mono: it ships on every Windows
    // back to Vista, and the monospace face carries no Fluent identity the way
    // Segoe UI Variable does — reliability is worth more here.
    readonly property string fontMono: "Consolas"
    // Segoe Fluent Icons is the Windows 11 system icon font. All codepoints
    // used here are shared with Segoe MDL2 Assets (Windows 10), so the glyphs
    // survive on older Windows even though the family name differs.
    readonly property string fontIcons: "Segoe Fluent Icons"

    // Canvas2D's ctx.font takes a CSS font shorthand string, which requires quoting
    // multi-word family names ("Segoe UI") or the parser silently mangles them into
    // "SegoeUI" (no such font) and falls back — QML's own font.family property has no
    // such requirement, so these quoted variants exist only for ctx.font assignments,
    // never for font.family: bindings.
    readonly property string fontFamilyCss: "'" + fontFamily + "'"
    readonly property string fontDisplayCss: "'" + fontDisplay + "'"
    readonly property string fontMonoCss: "'" + fontMono + "'"

    // Windows 11 type ramp. Fluent's Caption is 12px and Body is 14px; the old
    // scale ran 10/11/12, below the 12px floor for readable body text, which is
    // what made the ruler numerals and property labels feel cramped. Compact
    // density drops one step but never below Caption's 12px.
    readonly property real fontSizeCaption: 12   // Fluent Caption   12/16
    readonly property real fontSizeLabel: compactDensity ? 12 : 13
    readonly property real fontSizeBody: compactDensity ? 13 : 14   // Fluent Body 14/20
    readonly property real fontSizeHeadline: 16  // Fluent Body Large 18/24, tightened
    readonly property real fontSizeDisplay: 20   // Fluent Subtitle  20/28

    // Selection and Overlays
    readonly property color selectionOverlay: "#330078D4"
    readonly property color selectionOverlayStrong: "#660078D4"
    readonly property color dragRect: "#4FC3F7"
    readonly property color dragFill: "#0F4FC3F7"
    
    // ── Layout — Fluent 4px grid ───────────────────────────────────────────
    // Fluent's standard control height is 32px; its compact variant is 24px.
    // Density maps onto those two rather than inventing intermediate sizes.
    readonly property real controlHeight: compactDensity ? 24 : 32
    readonly property real iconSize: 16          // Fluent standard icon size
    readonly property real captionButtonWidth: 46   // Windows caption button
    readonly property real captionBarHeight: 32
    readonly property real toolbarHeight: compactDensity ? 32 : 40
    readonly property real menuBarHeight: 32
    // Tool button density (left panel, bottom atom bar, toolbar dropdowns)
    readonly property real toolCellSize: compactDensity ? 24 : 32
    readonly property real toolIconSize: 16
    readonly property real toolGridGap: 4
    // Wide enough for three toolCellSize columns + gaps + scrollbar
    readonly property real toolPanelWidth: compactDensity ? 96 : 120
    readonly property real propertyPanelWidth: compactDensity ? 208 : 248

    // Fluent corner radii: 4px on controls, 8px on cards/dialogs/flyouts.
    readonly property real radiusControl: 4
    readonly property real radiusSurface: 8

    // 4px spacing grid
    readonly property real marginSmall: 4
    readonly property real marginMedium: compactDensity ? 8 : 12
    readonly property real marginLarge: compactDensity ? 12 : 16
    readonly property real spacingSmall: 4
    readonly property real spacingMedium: 8
    readonly property real spacingLarge: 12

    // Dialog width scale — replaces the ad-hoc 300/360/380/420/480 spread.
    readonly property real dialogWidthSmall: 320   // single field + buttons
    readonly property real dialogWidthMedium: 420  // multi-field form
    readonly property real dialogWidthLarge: 520   // results, lists, tables
    readonly property real dialogHeightLarge: 560
    readonly property real dialogPadding: marginLarge
    readonly property real dialogSpacing: spacingLarge

    // ── Interaction states ──────────────────────────────────────────────────
    // Standard triad for icon-button opacity across enabled/hover/rest states.
    // Any new control should reference these rather than a literal — a few
    // older controls (e.g. MainWindow's transform/rotate-flip row) still use
    // one-off 0.7/0.25 values and should be reconciled to this triad when touched.
    readonly property color focusRingColor: accent
    readonly property real focusRingWidth: 2
    readonly property real disabledOpacity: 0.3
    readonly property real restOpacity: 0.6
    readonly property real hoverOpacity: 0.8

    // ── Motion (Fluent 2) ────────────────────────────────────────────────────
    // Curves sourced from Microsoft's own published Fluent motion guidance
    // (learn.microsoft.com/en-us/windows/apps/design/motion/timing-and-easing):
    // "Fast Out, Slow In" for anything entering, "Slow Out, Fast In" for anything
    // exiting. QML's Easing.BezierSpline takes a flat [c1x,c1y,c2x,c2y,endx,endy]
    // array -- the direct encoding of a CSS cubic-bezier(x1,y1,x2,y2) ending at (1,1).
    readonly property var easingDecelerate: [0, 0, 0, 1, 1, 1]   // cubic-bezier(0,0,0,1)
    readonly property var easingAccelerate: [1, 0, 1, 1, 1, 1]   // cubic-bezier(1,0,1,1)
    readonly property int durationFast: 150
    readonly property int durationMedium: 250
    readonly property int durationSlow: 400

    // Chemistry UI Constants (style-sheet driven) — proxied from StyleSheets.currentSheet.
    // Theme.qml owns UI chrome tokens (colors/typography/layout/interaction states above);
    // StyleSheets.qml owns user-selectable chemistry-rendering presets (ACD/ChemSketch, RSC, etc).
    // Don't add new literal overrides here — add new presets to StyleSheets.qml instead.
    property real atomRadius: StyleSheets.currentSheet.atomRadius
    property real bondWidth: StyleSheets.currentSheet.bondWidth
    property real selectionWidth: StyleSheets.currentSheet.selectionWidth
    property real hoverWidth: StyleSheets.currentSheet.hoverWidth
    property real baseFontSize: StyleSheets.currentSheet.baseFontSize
    property real baseBondLength: StyleSheets.currentSheet.baseBondLength

    // Stroke-width clamping (Phase 0): prevents lines from vanishing at low zoom
    // or becoming bloated at high zoom. Floor lowered from 0.8 to 0.4 (2026-07-18): the old
    // value happened to sit almost exactly at the old (unverified) ACS 1996 preset's line
    // width and would clamp RSC's real, verified 0.53px line thicker than the actual RSC
    // ChemDraw template specifies.
    readonly property real strokeWidthMin: 0.4
    readonly property real strokeWidthMax: 4.0
    property real wedgeTaperRatio: StyleSheets.currentSheet.wedgeTaperRatio
    property real doubleBondSpacing: StyleSheets.currentSheet.doubleBondSpacing
    readonly property real hashMinSize: 5
    property real hashSpacingFactor: StyleSheets.currentSheet.hashSpacingFactor
    property real wavyAmplitude: StyleSheets.currentSheet.wavyAmplitude
    readonly property real wavyWaveSpacing: 5.0
    readonly property real arrowHeadLength: 12.0
    readonly property real arrowHeadWidth: 5.0
    readonly property real plusSize: 6.0
    readonly property real subscriptOffset: 0.25
    readonly property real superscriptOffset: 0.3
    readonly property real fontSizeMin: 8.0
    readonly property real fontSizeMax: 24.0

    function clampStrokeWidth(width, scaleFactor) {
        return Math.max(strokeWidthMin, Math.min(strokeWidthMax, width * scaleFactor))
    }

    function clampDim(dim, scaleFactor) {
        return Math.max(0.5, dim * scaleFactor)
    }

    function clampFontSize(size, scaleFactor) {
        return Math.max(fontSizeMin, Math.min(fontSizeMax, size * scaleFactor))
    }

    readonly property var elementColors: {
        "nonmetal": "#A0E6FF",
        "alkali": "#FFB5B5",
        "alkaline": "#FFDCA8",
        "transition": "#FFC0CB",
        "post-transition": "#CCCCCC",
        "metalloid": "#99CC99",
        "halogen": "#FFFF99",
        "noble": "#FFE066",
        "lanthanide": "#FFBFFF",
        "actinide": "#FF99CC"
    }

    function getColorForType(type) {
        return elementColors[type] || background
    }

    readonly property var elementColorMap: { "H": "#000000", "C": "#000000", "N": "#3050F8", "O": "#FF0D0D", "F": "#90E050", "P": "#FF8000", "S": "#C39A00", "Cl": "#1FF01F", "Br": "#A62929", "I": "#940094" }

    function getElementColor(label) {
        if (!cpkColorsEnabled) return textPrimary;
        return elementColorMap[label] || textPrimary;
    }

    readonly property var elementsList: [
        {number:1,label:"H",type:"nonmetal"},{number:2,label:"He",type:"noble"},
        {number:3,label:"Li",type:"alkali"},{number:4,label:"Be",type:"alkaline"},
        {number:5,label:"B",type:"metalloid"},{number:6,label:"C",type:"nonmetal"},
        {number:7,label:"N",type:"nonmetal"},{number:8,label:"O",type:"nonmetal"},
        {number:9,label:"F",type:"halogen"},{number:10,label:"Ne",type:"noble"},
        {number:11,label:"Na",type:"alkali"},{number:12,label:"Mg",type:"alkaline"},
        {number:13,label:"Al",type:"post-transition"},{number:14,label:"Si",type:"metalloid"},
        {number:15,label:"P",type:"nonmetal"},{number:16,label:"S",type:"nonmetal"},
        {number:17,label:"Cl",type:"halogen"},{number:18,label:"Ar",type:"noble"},
        {number:19,label:"K",type:"alkali"},{number:20,label:"Ca",type:"alkaline"},
        {number:21,label:"Sc",type:"transition"},{number:22,label:"Ti",type:"transition"},
        {number:23,label:"V",type:"transition"},{number:24,label:"Cr",type:"transition"},
        {number:25,label:"Mn",type:"transition"},{number:26,label:"Fe",type:"transition"},
        {number:27,label:"Co",type:"transition"},{number:28,label:"Ni",type:"transition"},
        {number:29,label:"Cu",type:"transition"},{number:30,label:"Zn",type:"transition"},
        {number:31,label:"Ga",type:"post-transition"},{number:32,label:"Ge",type:"metalloid"},
        {number:33,label:"As",type:"metalloid"},{number:34,label:"Se",type:"nonmetal"},
        {number:35,label:"Br",type:"halogen"},{number:36,label:"Kr",type:"noble"},
        {number:37,label:"Rb",type:"alkali"},{number:38,label:"Sr",type:"alkaline"},
        {number:39,label:"Y",type:"transition"},{number:40,label:"Zr",type:"transition"},
        {number:41,label:"Nb",type:"transition"},{number:42,label:"Mo",type:"transition"},
        {number:43,label:"Tc",type:"transition"},{number:44,label:"Ru",type:"transition"},
        {number:45,label:"Rh",type:"transition"},{number:46,label:"Pd",type:"transition"},
        {number:47,label:"Ag",type:"transition"},{number:48,label:"Cd",type:"transition"},
        {number:49,label:"In",type:"post-transition"},{number:50,label:"Sn",type:"post-transition"},
        {number:51,label:"Sb",type:"metalloid"},{number:52,label:"Te",type:"metalloid"},
        {number:53,label:"I",type:"halogen"},{number:54,label:"Xe",type:"noble"},
        {number:55,label:"Cs",type:"alkali"},{number:56,label:"Ba",type:"alkaline"},
        {number:57,label:"La",type:"lanthanide"},{number:58,label:"Ce",type:"lanthanide"},
        {number:59,label:"Pr",type:"lanthanide"},{number:60,label:"Nd",type:"lanthanide"},
        {number:61,label:"Pm",type:"lanthanide"},{number:62,label:"Sm",type:"lanthanide"},
        {number:63,label:"Eu",type:"lanthanide"},{number:64,label:"Gd",type:"lanthanide"},
        {number:65,label:"Tb",type:"lanthanide"},{number:66,label:"Dy",type:"lanthanide"},
        {number:67,label:"Ho",type:"lanthanide"},{number:68,label:"Er",type:"lanthanide"},
        {number:69,label:"Tm",type:"lanthanide"},{number:70,label:"Yb",type:"lanthanide"},
        {number:71,label:"Lu",type:"lanthanide"},{number:72,label:"Hf",type:"transition"},
        {number:73,label:"Ta",type:"transition"},{number:74,label:"W",type:"transition"},
        {number:75,label:"Re",type:"transition"},{number:76,label:"Os",type:"transition"},
        {number:77,label:"Ir",type:"transition"},{number:78,label:"Pt",type:"transition"},
        {number:79,label:"Au",type:"transition"},{number:80,label:"Hg",type:"transition"},
        {number:81,label:"Tl",type:"post-transition"},{number:82,label:"Pb",type:"post-transition"},
        {number:83,label:"Bi",type:"post-transition"},{number:84,label:"Po",type:"post-transition"},
        {number:85,label:"At",type:"halogen"},{number:86,label:"Rn",type:"noble"},
        {number:87,label:"Fr",type:"alkali"},{number:88,label:"Ra",type:"alkaline"},
        {number:89,label:"Ac",type:"actinide"},{number:90,label:"Th",type:"actinide"},
        {number:91,label:"Pa",type:"actinide"},{number:92,label:"U",type:"actinide"},
        {number:93,label:"Np",type:"actinide"},{number:94,label:"Pu",type:"actinide"},
        {number:95,label:"Am",type:"actinide"},{number:96,label:"Cm",type:"actinide"},
        {number:97,label:"Bk",type:"actinide"},{number:98,label:"Cf",type:"actinide"},
        {number:99,label:"Es",type:"actinide"},{number:100,label:"Fm",type:"actinide"},
        {number:101,label:"Md",type:"actinide"},{number:102,label:"No",type:"actinide"},
        {number:103,label:"Lr",type:"actinide"},{number:104,label:"Rf",type:"transition"},
        {number:105,label:"Db",type:"transition"},{number:106,label:"Sg",type:"transition"},
        {number:107,label:"Bh",type:"transition"},{number:108,label:"Hs",type:"transition"},
        {number:109,label:"Mt",type:"transition"},{number:110,label:"Ds",type:"transition"},
        {number:111,label:"Rg",type:"transition"},{number:112,label:"Cn",type:"transition"},
        {number:113,label:"Nh",type:"post-transition"},{number:114,label:"Fl",type:"post-transition"},
        {number:115,label:"Mc",type:"post-transition"},{number:116,label:"Lv",type:"post-transition"},
        {number:117,label:"Ts",type:"halogen"},{number:118,label:"Og",type:"noble"}
    ]
}
