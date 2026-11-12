import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Rectangle {
    id: root
    property Item canvas
    // Set true only on the instance living inside floatingPropertyWindow
    // (MainWindow.qml). Purely cosmetic -- swaps which icon/tooltip/target-
    // state this panel's float/dock button shows; both instances otherwise
    // behave identically.
    property bool isFloating: false
    // Null-safe accessors: `canvas` is bound to activeCanvas, which is null for a
    // moment at startup before the per-document Repeater instantiates. Every binding
    // below reads these instead of dereferencing canvas directly.
    readonly property var selAtom: canvas ? canvas.selectedAtom : null
    readonly property var selBond: canvas ? canvas.selectedBond : null
    readonly property var selArrow: canvas ? canvas.selectedRxnArrow : null

    readonly property bool hasMultiSelection: {
        if (!canvas || !canvas.sketch) return false
        let atomCount = canvas.sketch.selection.atom_ids ? canvas.sketch.selection.atom_ids.length : 0
        let bondCount = canvas.sketch.selection.bond_ids ? canvas.sketch.selection.bond_ids.length : 0
        return (atomCount + bondCount) > 1
    }
    readonly property bool hasContextSelection: selAtom !== null || selBond !== null || selArrow !== null || hasMultiSelection
    readonly property var matchedRing: {
        if (!canvas || !canvas.sketch || !canvas.sketch.primitives || !canvas.sketch.primitives.rings) return null
        if (!canvas.sketch.selection || !canvas.sketch.selection.atom_ids) return null
        let selAtoms = canvas.sketch.selection.atom_ids
        let selBonds = canvas.sketch.selection.bond_ids || []
        if (selAtoms.length < 3) return null
        let selSet = {}
        for (let i = 0; i < selAtoms.length; i++) selSet[selAtoms[i]] = true
        let rings = canvas.sketch.primitives.rings
        for (let r = 0; r < rings.length; r++) {
            let ring = rings[r]
            if (!ring.atoms || ring.atoms.length !== selAtoms.length) continue
            if (selBonds.length !== ring.atoms.length) continue
            let match = true
            for (let j = 0; j < ring.atoms.length; j++) {
                if (!selSet[ring.atoms[j]]) { match = false; break }
            }
            if (match) return ring
        }
        return null
    }
    property bool drugPropertiesExpanded: true

    property string molName: ""
    property var sdfProps: ({})
    property double molMW:       0
    property double molMono:     0
    property string molFormula:  ""
    property int    molAtoms:    0
    property int    molBonds:    0
    property double molTPSA:     0
    property double molLogP:     0
    property double molMolarRefractivity: 0
    property double molPka: 0
    property int    molHeavyAtoms: 0
    property bool   molIsChiral: false
    property double molMostAbundantMass: 0
    property int    molFragmentCount: 0
    property int    molRingCount: 0
    property int    molHBA:      0
    property int    molHBD:      0
    property int    molRotBonds: 0

    property string iupacName:    ""
    property string iupacError:   ""
    property bool   iupacLoading: false

    // This panel is a single window-global instance shared across every open
    // document -- without this, switching tabs left the previous document's
    // generated name showing under the newly-active (different) molecule.
    onCanvasChanged: {
        iupacName = ""
        iupacError = ""
        iupacLoading = false
    }

    // Indigo's C API returns exactly -1 as a sentinel when a scalar property
    // calculation fails internally. Failures can be PARTIAL within one
    // calcProperties batch -- observed live with an invalid-valence structure:
    // MW came back 196.21 (computed from the formula) while Exact mass, TPSA,
    // MolRef, pKa and HBD all returned -1 simultaneously. So a whole-batch
    // check on MW alone misses real failures. Exact (monoisotopic) mass can
    // never legitimately be negative, so molMW/molMono double as the batch
    // signal; every field that cannot be negative is additionally guarded
    // per-value. LogP and pKa CAN be legitimately negative (caffeine LogP is
    // -1.03), so they rely on the batch signal only.
    readonly property bool calcFailed: root.molMW < 0 || root.molMono < 0

    function _fmtNum(v, decimals, suffix, canBeNegative) {
        if (root.calcFailed || (!canBeNegative && v < 0)) return "—"
        return v.toFixed(decimals) + (suffix || "")
    }
    function _fmtInt(v) { return (root.calcFailed || v < 0) ? "—" : v.toString() }
    function _formulaHtml(f) { return f.replace(/([0-9]+)/g, "<sub>$1</sub>") }
    function _panelTitle() {
        if (root.selAtom !== null) return "Atom \u00b7 " + (root.selAtom.label || "?") + root.selAtom.id
        if (root.selBond !== null) return "Bond \u00b7 #" + root.selBond.id
        if (root.selArrow !== null) return "Reaction"
        if (root.matchedRing !== null) return "Ring \u00b7 " + root.matchedRing.atoms.length + "-membered"
        if (root.hasMultiSelection) {
            let atomCount = (root.canvas && root.canvas.sketch && root.canvas.sketch.selection.atom_ids) ? root.canvas.sketch.selection.atom_ids.length : 0
            let bondCount = (root.canvas && root.canvas.sketch && root.canvas.sketch.selection.bond_ids) ? root.canvas.sketch.selection.bond_ids.length : 0
            return "Selection \u00b7 " + (atomCount + bondCount) + " objects"
        }
        if (root.molAtoms > 0) return "Molecule"
        return "Properties"
    }

    // Indigo reports check types as raw keys ("ambiguous_h", "overlap_atom").
    // Shown verbatim they read as defects in the drawing rather than what they
    // are, which cost the whole checker credibility -- "ambiguous h" on a
    // correctly drawn caffeine is Indigo saying the aromatic ring form leaves
    // the implicit-H count undetermined, not that the structure is wrong.
    readonly property var _issueLabels: ({
        "ambiguous_h":  "implicit H count undetermined (aromatic form)",
        "valence":      "unusual valence",
        "radical":      "radical centre",
        "pseudoatom":   "pseudoatom",
        "stereo":       "stereo descriptor problem",
        "3d_coord":     "has 3D coordinates",
        "overlap_atom": "overlapping atoms",
        "overlap_bond": "overlapping bonds",
        "query":        "query feature",
        "chirality":    "chiral flag inconsistent",
        "chiral_flag":  "chiral flag inconsistent"
    })
    function _issueText(w) { return root._issueLabels[w] || w.replace(/_/g, " ") }

    color: Theme.surface

    ColumnLayout {
        anchors {
            fill: parent
            margins: Theme.marginLarge
        }
        spacing: Theme.spacingLarge

        RowLayout {
            Layout.fillWidth: true
            spacing: Theme.spacingSmall

            Text {
                textFormat: Text.PlainText
                text: root._panelTitle()
                color: Theme.textSecondary
                font {
                    pixelSize: Theme.fontSizeCaption
                    bold: true
                    letterSpacing: 2
                    family: Theme.fontFamily
                }
                Layout.fillWidth: true
                horizontalAlignment: Text.AlignHCenter
            }

            // Collapse: reclaims this panel's full width back to the canvas.
            // Re-expand happens via the slim tab MainWindow shows in its place
            // (Step 3) or the bottom bar's issue badge (Task 2).
            ToolButton {
                // Theme.fontIcons (Segoe Fluent Icons) private-use-area glyphs
                // do not reliably survive every editing path -- confirmed by
                // screenshot that this exact pattern renders empty elsewhere
                // in this codebase (MainWindow.qml's density toggle, same
                // issue, pre-existing). Plain font + a standard BMP arrow
                // avoids the whole class of bug; matches the re-expand tab's
                // "◀" below, which is already proven to render correctly.
                text: "▶"
                font { family: Theme.fontFamily; pixelSize: Theme.iconSize }
                implicitWidth: Theme.controlHeight
                implicitHeight: Theme.controlHeight
                Accessible.name: "Collapse properties panel"
                ToolTip.text: "Collapse"
                ToolTip.visible: hovered
                ToolTip.delay: 500
                onClicked: uiSettings.panelState = "collapsed"
                visible: !root.isFloating
            }
            // Float: detaches this panel's content into its own draggable
            // window (floatingPropertyWindow in MainWindow.qml). The SAME
            // button, shown in BOTH the docked and floating instances, docks
            // back when clicked while already floating -- root.isFloating
            // (set by MainWindow when instantiating the floating instance)
            // tells each instance which state to toggle TO on click.
            ToolButton {
                // Same private-use-glyph issue as the collapse button above --
                // plain font + standard arrows instead.
                text: root.isFloating ? "↙" : "↗"
                font { family: Theme.fontFamily; pixelSize: Theme.iconSize }
                implicitWidth: Theme.controlHeight
                implicitHeight: Theme.controlHeight
                Accessible.name: root.isFloating ? "Dock properties panel" : "Float properties panel"
                ToolTip.text: root.isFloating ? "Dock" : "Float"
                ToolTip.visible: hovered
                ToolTip.delay: 500
                onClicked: uiSettings.panelState = root.isFloating ? "docked" : "floating"
            }
        }

        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 1
            color: Theme.outline
        }

        // Molecular Properties (auto-updated)
        ColumnLayout {
            Layout.fillWidth: true
            visible: root.molAtoms > 0 && !root.hasContextSelection
            spacing: 6

            TextField {
                id: molNameInput
                Layout.fillWidth: true
                placeholderText: "Untitled"
                text: root.molName
                font { pixelSize: Theme.fontSizeBody; family: Theme.fontFamily }
                onEditingFinished: {
                    if (root.canvas && root.canvas.sketch) {
                        root.canvas.sketch.sendCommand("setMoleculeName", [text])
                    }
                    focus = false
                }
                Keys.onEscapePressed: {
                    text = root.molName
                    focus = false
                }
            }

            // IUPAC Name (beta) section. Label and button stack on separate
            // rows rather than sharing one -- Fluent's Button carries wider
            // horizontal padding than the old Basic style, and at this panel's
            // width that left almost nothing for the label, which is why it
            // was truncating to "IUPAC Na..." even with elide already in place.
            ColumnLayout {
                Layout.fillWidth: true
                spacing: 4

                Text {
                    textFormat: Text.PlainText
                    text: "IUPAC Name (beta)"
                    color: Theme.textSecondary
                    font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
                    Layout.fillWidth: true
                }

                Button {
                    Layout.alignment: Qt.AlignRight
                    text: root.iupacLoading ? "Generating…" : "Generate"
                    // ChemDraw's "Structure to Name" convention: selection required, named
                    // only when something is selected (never silently names the whole
                    // document -- a multi-structure canvas would otherwise be ambiguous
                    // about which structure the name even describes).
                    enabled: !root.iupacLoading && !!(root.canvas && root.canvas.sketch &&
                             root.canvas.sketch.selection && root.canvas.sketch.selection.atom_ids &&
                             root.canvas.sketch.selection.atom_ids.length > 0)
                    font { pixelSize: Theme.fontSizeCaption; family: Theme.fontFamily }
                    ToolTip.visible: hovered && !enabled && !root.iupacLoading
                    ToolTip.text: "Select the structure to name first"
                    onClicked: {
                        if (root.canvas && root.canvas.sketch) {
                            root.iupacLoading = true
                            root.iupacName = ""
                            root.iupacError = ""
                            root.canvas.sketch.requestMolfileForNamingSelection("iupac_name")
                        }
                    }
                }

                Rectangle {
                    Layout.fillWidth: true
                    implicitHeight: iupacText.implicitHeight + 10
                    visible: root.iupacName !== "" || root.iupacError !== ""
                    color: root.iupacError !== "" ? "#22FF0000" : Theme.surface
                    border.color: root.iupacError !== "" ? "#66FF0000" : Theme.outline
                    radius: 4

                    Text {
                        id: iupacText
                        anchors { fill: parent; margins: 5 }
                        wrapMode: Text.Wrap
                        textFormat: Text.PlainText
                        text: root.iupacError !== "" ? root.iupacError : root.iupacName
                        color: root.iupacError !== "" ? "#FF6B6B" : Theme.textPrimary
                        font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono; bold: root.iupacName !== "" }
                    }
                }
            }

            Text {
                textFormat: Text.RichText
                text: root._formulaHtml(root.molFormula)
                color: Theme.textPrimary
                // Deliberately the one serif hold-out in this panel: the molecular
                // formula is the chemistry, and reads as a journal figure caption.
                // Surrounding chrome labels are sans (Theme.fontFamily).
                font {
                    pixelSize: Theme.fontSizeDisplay
                    bold: true
                    family: Theme.fontDisplay
                }
                Layout.alignment: Qt.AlignHCenter
            }

            Grid {
                columns: 2
                columnSpacing: 12
                rowSpacing: 4
                Layout.fillWidth: true

                Text { textFormat: Text.PlainText; text: "MW";        color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtNum(root.molMW, 3, " g/mol"); color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }
                Text { textFormat: Text.PlainText; text: "Exact";     color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtNum(root.molMono, 4, " Da");   color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }
                Text { textFormat: Text.PlainText; text: "Atoms";     color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root.molAtoms.toString();                 color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }
                Text { textFormat: Text.PlainText; text: "Bonds";     color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root.molBonds.toString();                 color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }
            }

            Rectangle { Layout.fillWidth: true; height: 1; color: Theme.outline; opacity: 0.5 }

            // MouseArea kept outside the RowLayout: anchors don't apply to layout children.
            Item {
                Layout.fillWidth: true
                implicitHeight: drugHeaderRow.implicitHeight

                RowLayout {
                    id: drugHeaderRow
                    anchors.fill: parent
                    spacing: 4

                    Text {
                        textFormat: Text.PlainText
                        text: "Drug Properties"
                        color: Theme.textSecondary
                        font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
                        Layout.fillWidth: true
                    }
                    Text {
                        textFormat: Text.PlainText
                        text: root.drugPropertiesExpanded ? "\u25be" : "\u25b8"
                        color: Theme.textSecondary
                        font { pixelSize: Theme.fontSizeCaption; family: Theme.fontFamily }
                    }
                }

                MouseArea {
                    anchors.fill: parent
                    cursorShape: Qt.PointingHandCursor
                    onClicked: root.drugPropertiesExpanded = !root.drugPropertiesExpanded
                }
            }

            Grid {
                columns: 2
                columnSpacing: 12
                rowSpacing: 4
                Layout.fillWidth: true
                visible: root.drugPropertiesExpanded

                Text { textFormat: Text.PlainText; text: "TPSA";   color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtNum(root.molTPSA, 2, " Å²");  color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "LogP";   color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtNum(root.molLogP, 2, "", true);  color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "HBA";    color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtInt(root.molHBA);                color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "HBD";    color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtInt(root.molHBD);                color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "RotB";   color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtInt(root.molRotBonds);           color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "MolRef";  color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtNum(root.molMolarRefractivity, 2);  color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "pKa";     color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtNum(root.molPka, 2, "", true);  color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "Heavy Atoms"; color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtInt(root.molHeavyAtoms); color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "Chiral"; color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root.molIsChiral ? "Yes" : "No"; color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "Abundant Mass"; color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root._fmtNum(root.molMostAbundantMass, 3, " g/mol"); color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "Fragments"; color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root.molFragmentCount.toString(); color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "Rings (SSSR)"; color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root.molRingCount.toString(); color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }
            }
        }

        // SDF Data Fields (read-only)
        ColumnLayout {
            Layout.fillWidth: true
            visible: Object.keys(root.sdfProps).length > 0
            spacing: 6

            Rectangle { Layout.fillWidth: true; height: 1; color: Theme.outline; opacity: 0.5 }

            Text {
                textFormat: Text.PlainText
                text: "SDF Data Fields"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }

            Grid {
                columns: 2
                columnSpacing: 12
                rowSpacing: 4
                Layout.fillWidth: true

                Repeater {
                    model: {
                        var keys = Object.keys(root.sdfProps);
                        var arr = [];
                        for (var i = 0; i < keys.length; i++) {
                            arr.push({ isKey: true, val: keys[i] });
                            arr.push({ isKey: false, val: root.sdfProps[keys[i]] });
                        }
                        return arr;
                    }
                    delegate: Text {
                        textFormat: Text.PlainText
                        text: modelData.val
                        color: modelData.isKey ? Theme.textSecondary : Theme.textPrimary
                        font {
                            pixelSize: Theme.fontSizeLabel
                            family: modelData.isKey ? Theme.fontDisplay : Theme.fontMono
                        }
                    }
                }
            }
        }

        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 1
            color: Theme.outline
            visible: root.molAtoms > 0
        }

        // Multi-select summary
        ColumnLayout {
            Layout.fillWidth: true
            visible: root.hasMultiSelection && root.matchedRing === null
            spacing: Theme.spacingMedium

            Text {
                textFormat: Text.PlainText
                text: "Selection Summary"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }
            Text {
                textFormat: Text.PlainText
                text: {
                    if (!root.canvas || !root.canvas.sketch) return "";
                    let atomCount = root.canvas.sketch.selection.atom_ids ? root.canvas.sketch.selection.atom_ids.length : 0;
                    let bondCount = root.canvas.sketch.selection.bond_ids ? root.canvas.sketch.selection.bond_ids.length : 0;
                    let parts = [];
                    if (atomCount > 0) parts.push(atomCount + (atomCount === 1 ? " atom" : " atoms"));
                    if (bondCount > 0) parts.push(bondCount + (bondCount === 1 ? " bond" : " bonds"));
                    return parts.join(", ") + " selected";
                }
                color: Theme.textPrimary
                font { pixelSize: Theme.fontSizeBody; family: Theme.fontFamily }
            }
        }

        // Ring Properties
        ColumnLayout {
            Layout.fillWidth: true
            visible: root.matchedRing !== null
            spacing: Theme.spacingMedium

            Text {
                textFormat: Text.PlainText
                text: "Ring"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }

            Grid {
                columns: 2
                columnSpacing: 12
                rowSpacing: 4
                Layout.fillWidth: true

                Text { textFormat: Text.PlainText; text: "Size";     color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root.matchedRing ? (root.matchedRing.atoms.length + "-membered") : ""; color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "Aromatic"; color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: (root.matchedRing && root.matchedRing.isAromatic) ? "Yes" : "No"; color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }

                Text { textFormat: Text.PlainText; text: "Atoms";    color: Theme.textSecondary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily } }
                Text { textFormat: Text.PlainText; text: root.matchedRing ? root.matchedRing.atoms.length.toString() : ""; color: Theme.textPrimary; font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono } }
            }
        }

        // Check Issues Section
        ColumnLayout {
            Layout.fillWidth: true
            spacing: Theme.spacingMedium

            property var issuesList: {
                if (!root.canvas || !root.canvas.sketch || !root.canvas.sketch.primitives) return []
                var prims = root.canvas.sketch.primitives
                var lst = []
                if (prims.atoms) {
                    for (var i = 0; i < prims.atoms.length; i++) {
                        if (prims.atoms[i].checkWarning) {
                            lst.push({ type: "atom", id: prims.atoms[i].id, warning: prims.atoms[i].checkWarning })
                        }
                    }
                }
                if (prims.bonds) {
                    for (var j = 0; j < prims.bonds.length; j++) {
                        if (prims.bonds[j].checkWarning) {
                            lst.push({ type: "bond", id: prims.bonds[j].id, warning: prims.bonds[j].checkWarning })
                        }
                    }
                }
                return lst
            }

            visible: issuesList.length > 0

            Rectangle { Layout.fillWidth: true; height: 1; color: Theme.outline; opacity: 0.5 }

            RowLayout {
                Layout.fillWidth: true
                Text {
                    textFormat: Text.PlainText
                    text: "Check Issues"
                    color: Theme.textSecondary
                    font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
                    Layout.fillWidth: true
                }
                Text {
                    text: parent.parent.issuesList.length.toString()
                    color: Theme.error
                    font { pixelSize: Theme.fontSizeCaption; bold: true; family: Theme.fontFamily }
                }
            }

            ListView {
                Layout.fillWidth: true
                implicitHeight: Math.min(200, contentHeight)
                clip: true
                model: parent.issuesList
                spacing: 4
                delegate: Rectangle {
                    width: ListView.view.width
                    height: 24
                    color: mouseArea.containsMouse ? Theme.hover : "transparent"
                    radius: 4

                    Text {
                        anchors.verticalCenter: parent.verticalCenter
                        anchors.left: parent.left
                        anchors.leftMargin: 4
                        anchors.right: parent.right
                        text: (modelData.type === "atom" ? "Atom " : "Bond ") + modelData.id + ": " + root._issueText(modelData.warning)
                        color: Theme.error
                        font { pixelSize: Theme.fontSizeLabel; family: Theme.fontFamily }
                        elide: Text.ElideRight
                    }
                    MouseArea {
                        id: mouseArea
                        anchors.fill: parent
                        hoverEnabled: true
                        cursorShape: Qt.PointingHandCursor
                        onClicked: {
                            if (root.canvas && root.canvas.sketch) {
                                // set before selectItem(), which updates the canvas synchronously
                                root.canvas._needsCentering = true
                                if (modelData.type === "atom") {
                                    root.canvas.sketch.selectItem(modelData.id, null)
                                } else {
                                    root.canvas.sketch.selectItem(null, modelData.id)
                                }
                            }
                        }
                    }
                }
            }
        }

        // Sgroup (contracted abbreviation) info
        ColumnLayout {
            Layout.fillWidth: true
            visible: root.selAtom !== null && root.selAtom.isSgroup
            spacing: Theme.spacingMedium

            Text {
                textFormat: Text.PlainText
                text: "Abbreviation"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }
            TextField {
                id: sgroupLabelInput
                Layout.fillWidth: true
                text: root.selAtom ? root.selAtom.label : ""
                font { pixelSize: Theme.fontSizeHeadline; bold: true; family: Theme.fontFamily }
                onEditingFinished: {
                    if (root.selAtom && root.selAtom.isSgroup) {
                        canvas.applyPropertyChange("sgroupLabel", root.selAtom.id, text)
                    }
                    focus = false
                }
                Keys.onEscapePressed: {
                    text = root.selAtom ? root.selAtom.label : ""
                    focus = false
                }
            }
        }

        // Atom Properties (real atoms only)
        ColumnLayout {
            Layout.fillWidth: true
            visible: root.selAtom !== null && !root.selAtom.isSgroup
            spacing: Theme.spacingMedium

            Text {
                textFormat: Text.PlainText
                text: "Atom Label"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }
            TextField {
                id: atomLabelInput
                Layout.fillWidth: true
                text: root.selAtom ? root.selAtom.label : ""
                font { pixelSize: Theme.fontSizeBody; family: Theme.fontMono }
                onEditingFinished: {
                    if (root.selAtom && !root.selAtom.isSgroup) {
                        canvas.applyPropertyChange("atomLabel", root.selAtom.id, text)
                    }
                    focus = false
                }
                Keys.onEscapePressed: {
                    text = root.selAtom ? root.selAtom.label : ""
                    focus = false
                }
            }

            Text {
                textFormat: Text.PlainText
                text: "Charge"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }
            TextField {
                id: atomChargeInput
                Layout.fillWidth: true
                text: root.selAtom ? (root.selAtom.charge > 0 ? "+" + root.selAtom.charge : root.selAtom.charge.toString()) : "0"
                font { pixelSize: Theme.fontSizeBody; family: Theme.fontMono }
                onEditingFinished: {
                    if (root.selAtom && !root.selAtom.isSgroup) {
                        canvas.applyPropertyChange("atomCharge", root.selAtom.id, text)
                    }
                    focus = false
                }
                Keys.onEscapePressed: {
                    text = root.selAtom ? (root.selAtom.charge > 0 ? "+" + root.selAtom.charge : root.selAtom.charge.toString()) : "0"
                    focus = false
                }
            }

            Text {
                textFormat: Text.PlainText
                text: "Atom List (e.g. C,N,O)"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }
            RowLayout {
                Layout.fillWidth: true
                TextField {
                    id: atomListInput
                    Layout.fillWidth: true
                    placeholderText: "leave empty to clear"
                    text: root.selAtom ? root.selAtom.atomListElements : ""
                    font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono }
                    onEditingFinished: {
                        if (root.selAtom && !root.selAtom.isSgroup) {
                            canvas.applyPropertyChange("atomQueryList", root.selAtom.id,
                                { elements: text, notList: atomListNotCheck.checked })
                        }
                        focus = false
                    }
                    Keys.onEscapePressed: {
                        text = root.selAtom ? root.selAtom.atomListElements : ""
                        focus = false
                    }
                }
                CheckBox {
                    id: atomListNotCheck
                    text: "NOT"
                    checked: root.selAtom ? root.selAtom.atomListNot : false
                    onToggled: {
                        if (root.selAtom && !root.selAtom.isSgroup && atomListInput.text.trim().length > 0) {
                            canvas.applyPropertyChange("atomQueryList", root.selAtom.id,
                                { elements: atomListInput.text, notList: checked })
                        }
                    }
                }
            }
        }

        // Bond Properties
        ColumnLayout {
            Layout.fillWidth: true
            visible: root.selBond !== null
            spacing: Theme.spacingMedium

            Text {
                textFormat: Text.PlainText
                text: "Bond Type"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }
            ComboBox {
                id: bondTypeCombo
                Layout.fillWidth: true
                model: ["Single", "Double", "Triple"]
                onActivated: {
                    if (root.selBond) {
                        canvas.applyPropertyChange("bondType", root.selBond.id, currentIndex + 1)
                    }
                }
                Binding on currentIndex {
                    value: root.selBond ? Math.max(0, Math.min(2, root.selBond.type - 1)) : 0
                    restoreMode: Binding.RestoreBinding
                }
            }
        }

        // Reaction Arrow Properties
        ColumnLayout {
            Layout.fillWidth: true
            visible: root.selArrow !== null
            spacing: Theme.spacingMedium

            Text {
                textFormat: Text.PlainText
                text: "Arrow Mode"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }
            ComboBox {
                id: arrowModeCombo
                Layout.fillWidth: true
                model: ["Filled Triangle", "Open Angle", "Retrosynthetic", "Equilibrium (↔)", "Equilibrium (⇌)", "Curved Mechanism"]
                onActivated: {
                    if (root.selArrow) {
                        const modes = ["filled-triangle", "open-angle", "retrosynthetic", "equilibrium-FF", "equilibrium-FH", "curved-mechanism"]
                        canvas.applyPropertyChange("rxnArrowMode", root.selArrow.id, modes[currentIndex])
                    }
                }
                Binding on currentIndex {
                    value: {
                        if (!root.selArrow) return 0
                        const map = {"filled-triangle":0, "open-angle":1, "retrosynthetic":2, "equilibrium-FF":3, "equilibrium-FH":4, "curved-mechanism":5}
                        return map[root.selArrow.mode] !== undefined ? map[root.selArrow.mode] : 0
                    }
                    restoreMode: Binding.RestoreBinding
                }
            }

            Text {
                textFormat: Text.PlainText
                text: "Conditions Above"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }
            TextField {
                id: condAboveField
                Layout.fillWidth: true
                placeholderText: "e.g. 25°C, Pd/C"
                font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono }
                onEditingFinished: {
                    if (root.selArrow) {
                        canvas.applyPropertyChange("rxnArrowConditions", root.selArrow.id, {above: text, below: condBelowField.text})
                    }
                    focus = false
                }
                Keys.onEscapePressed: {
                    text = root.selArrow ? root.selArrow.conditionsAbove : ""
                    focus = false
                }
                Binding on text {
                    value: root.selArrow ? root.selArrow.conditionsAbove : ""
                    restoreMode: Binding.RestoreBinding
                }
            }

            Text {
                textFormat: Text.PlainText
                text: "Conditions Below"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontFamily }
            }
            TextField {
                id: condBelowField
                Layout.fillWidth: true
                placeholderText: "e.g. EtOH, 12h"
                font { pixelSize: Theme.fontSizeLabel; family: Theme.fontMono }
                onEditingFinished: {
                    if (root.selArrow) {
                        canvas.applyPropertyChange("rxnArrowConditions", root.selArrow.id, {above: condAboveField.text, below: text})
                    }
                    focus = false
                }
                Keys.onEscapePressed: {
                    text = root.selArrow ? root.selArrow.conditionsBelow : ""
                    focus = false
                }
                Binding on text {
                    value: root.selArrow ? root.selArrow.conditionsBelow : ""
                    restoreMode: Binding.RestoreBinding
                }
            }
        }

        // Empty state — centered in whatever space remains, matching the
        // canvas's own "Empty Canvas" placeholder convention. Previously this
        // was top-anchored followed by a bare filler Item, which left the rest
        // of the panel (up to ~700px when no molecule is loaded) as dead blank
        // space rather than an intentional empty state.
        Item {
            Layout.fillWidth: true
            Layout.fillHeight: true
            StatusPlaceholder {
                anchors.centerIn: parent
                visible: !root.hasContextSelection && root.molAtoms === 0
                mode: "empty"
                message: "Nothing selected"
                detail: "Click an atom, bond, or arrow\nto edit its properties"
            }
        }

        Text {
            textFormat: Text.PlainText
            text: "sketch"
            color: Theme.textSecondary
            opacity: 0.5
            font { pixelSize: Theme.fontSizeCaption; family: Theme.fontFamily; italic: true }
        }
    }
}
