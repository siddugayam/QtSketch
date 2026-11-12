import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import "js/Selection.js" as Selection

// Extracted from MainWindow.qml (was the inline AppMenuBar{...} instance). Takes
// the root ApplicationWindow via `win` (a required property, not an id lookup) —
// every reference to window state/functions/dialogs below is qualified through
// `win.` rather than bare, since bare names declared on `window` are not
// resolvable here across this new file boundary the same way they aren't
// resolvable across the AppMenuBarItem Popup-reparenting boundary that caused
// the earlier executeStructureOp ReferenceError bug (see window.executeStructureOp's
// own comment in MainWindow.qml for the full explanation of why qmlsc AOT needs
// qualified paths here).
AppMenuBar {
    id: root
    required property var win

    AppMenuBarItem {
        text: "File"
        MenuItemRow { text: "New Document"; onTriggered: DocumentManager.addDocument() }
        MenuItemRow { text: "Open…"; iconSource: "open.svg"; onTriggered: root.win.openDialog.open() }
        MenuItemRow { text: "Save"; iconSource: "save.svg"; shortcutHint: "Ctrl+S"; onTriggered: root.win.saveActive(false) }
        MenuItemRow { text: "Save As…"; iconSource: "save.svg"; onTriggered: root.win.saveActive(true) }

        // Recent Files
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6; visible: root.win.uiSettings.recentFilesJoined !== "" }
        Repeater {
            model: root.win.uiSettings.recentFilesJoined ? root.win.uiSettings.recentFilesJoined.split("|") : []
            MenuItemRow {
                text: modelData.substring(modelData.lastIndexOf("/") + 1) || modelData
                iconSource: "open.svg"
                onTriggered: root.win.loadFromFile(Qt.url(modelData))
            }
        }

        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Load from SMILES…"; iconSource: "smiles_in.svg"; onTriggered: root.win.smilesDialog.open() }
        MenuItemRow { text: "Load from InChI…"; iconSource: "smiles_in.svg"; onTriggered: root.win.inchiLoadDialog.open() }
        MenuItemRow { text: "Biopolymer…"; iconSource: "biopolymer.svg"; onTriggered: root.win.biopolymerDialog.open() }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Export as PDF…"; iconSource: "file-thumbnail.svg"; onTriggered: root.win.printToPdf() }
        MenuItemRow { text: "Export Reaction Scheme (Grid)…"; iconSource: "file-thumbnail.svg"; onTriggered: root.win.gridSaveDialog.open() }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Clear Canvas"; iconSource: "clear.svg"; onTriggered: root.win.activeCanvas.clearCanvas() }
    }
    AppMenuBarItem {
        text: "Edit"
        MenuItemRow { text: "Undo"; iconSource: "undo.svg"; shortcutHint: "Ctrl+Z"; enabled: !!(root.win.activeCanvas && root.win.activeCanvas.canUndo); onTriggered: root.win.activeCanvas.undo() }
        MenuItemRow { text: "Redo"; iconSource: "redo.svg"; shortcutHint: "Ctrl+Y"; enabled: !!(root.win.activeCanvas && root.win.activeCanvas.canRedo); onTriggered: root.win.activeCanvas.redo() }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Cut"; iconSource: "cut.svg"; shortcutHint: "Ctrl+X"; onTriggered: root.win.activeCanvas.cutSelection() }
        MenuItemRow { text: "Copy"; iconSource: "copy.svg"; shortcutHint: "Ctrl+C"; onTriggered: root.win.activeCanvas.copySelection() }
        MenuItemRow { text: "Paste"; iconSource: "paste.svg"; shortcutHint: "Ctrl+V"; onTriggered: root.win.activeCanvas.pasteSelection() }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Copy as Image"; iconSource: "copy_image.svg"; shortcutHint: "Ctrl+Shift+C"; onTriggered: root.win.activeCanvas.copyAsImage() }
        MenuItemRow { text: "Save Selection as Template…"; iconSource: "copy.svg"; enabled: Selection.hasAtoms(root.win.activeSketch, 1); onTriggered: root.win.saveSelectionAsTemplate() }
    }
    AppMenuBarItem {
        text: "Structure"
        MenuItemRow { text: "Layout"; iconSource: "layout.svg"; onTriggered: root.win.executeStructureOp("layout") }
        MenuItemRow { text: "Layout Selected"; iconSource: "layout.svg"; enabled: Selection.hasAtoms(root.win.activeSketch, 1) && !root.win.isProcessing; onTriggered: root.win.executeStructureOp("layoutSelected") }
        MenuItemRow { text: "Clean 2D"; iconSource: "layout.svg"; onTriggered: root.win.executeStructureOp("clean2d") }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Aromatize"; iconSource: "arom.svg"; onTriggered: root.win.executeStructureOp("aromatize") }
        MenuItemRow { text: "Dearomatize"; iconSource: "dearom.svg"; onTriggered: root.win.executeStructureOp("dearomatize") }
        MenuItemRow { text: "Add Explicit H"; iconSource: "explicit-hydrogens.svg"; onTriggered: root.win.executeStructureOp("unfoldH") }
        MenuItemRow { text: "Remove Explicit H"; iconSource: "explicit-hydrogens.svg"; onTriggered: root.win.executeStructureOp("foldH") }
        MenuItemRow { text: "Normalize"; iconSource: "clean.svg"; onTriggered: root.win.executeStructureOp("normalize") }
        MenuItemRow { text: "Standardize"; iconSource: "analyse.svg"; onTriggered: root.win.executeStructureOp("standardize") }
        MenuItemRow { text: "Ionize at pH…"; iconSource: "analyse.svg"; onTriggered: root.win.ionizeDialog.open() }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Validate"; iconSource: "check.svg"; onTriggered: root.win.executeStructureOp("check") }
        MenuItemRow { text: "Search Substructure (SMARTS)…"; iconSource: "search.svg"; onTriggered: root.win.executeStructureOp("smartsSearch") }
        MenuItemRow { text: "R-Groups…"; iconSource: "rgroup-label.svg"; shortcutHint: "Ctrl+R"; onTriggered: root.win.executeStructureOp("rgroups") }
    }
    AppMenuBarItem {
        id: reactionsMenuTrigger
        text: "Reactions"
        ColumnLayout {
            width: 210
            spacing: 6

            GridLayout {
                columns: 4
                columnSpacing: Theme.spacingSmall
                Layout.alignment: Qt.AlignHCenter
                Repeater {
                    model: [
                        { id: "RXN_ARROW", icon: "reaction-arrow-open-angle.svg", tip: "Reaction Arrow" },
                        { id: "MULTITAIL_ARROW", icon: "reaction-arrow-multitail.svg", tip: "Multi-tail Arrow" },
                        { id: "RXN_PLUS", icon: "reaction-plus.svg", tip: "Reaction Plus" },
                        { id: "AAM", icon: "reaction-map.svg", tip: "Atom-Atom Mapping" }
                    ]
                    delegate: IconCell {
                        id: rxnToolCell
                        required property var modelData
                        cellWidth: 44
                        cellHeight: Theme.toolCellSize
                        iconSource: "icons/" + rxnToolCell.modelData.icon
                        tip: rxnToolCell.modelData.tip
                        selected: !!(root.win.activeCanvas && root.win.activeCanvas.currentTool === rxnToolCell.modelData.id)
                        onClicked: {
                            if (!root.win.activeCanvas) return
                            root.win.activeCanvas.currentTool = root.win.activeCanvas.currentTool === rxnToolCell.modelData.id ? "SELECT" : rxnToolCell.modelData.id
                            reactionsMenuTrigger.close()
                        }
                    }
                }
            }

            Rectangle { Layout.fillWidth: true; Layout.preferredHeight: 1; color: Theme.outline }

            Text {
                text: "Arrow style"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontDisplay }
            }
            Repeater {
                model: [
                    { mode: "filled-triangle", label: "Filled Triangle", icon: "reaction-arrow-filled-triangle.svg" },
                    { mode: "open-angle", label: "Open Angle", icon: "reaction-arrow-open-angle.svg" },
                    { mode: "retrosynthetic", label: "Retrosynthetic", icon: "reaction-arrow-retrosynthetic-arrow.svg" },
                    { mode: "equilibrium-FF", label: "Equilibrium (↔)", icon: "reaction-arrow-equilibrium-filled-triangle.svg" },
                    { mode: "equilibrium-FH", label: "Equilibrium (⇌)", icon: "reaction-arrow-equilibrium-filled-half-bow.svg" },
                    { mode: "curved-mechanism", label: "Curved Mechanism", icon: "reaction-arrow-elliptical-arc-arrow-filled-triangle.svg" }
                ]
                delegate: Rectangle {
                    id: amRow
                    required property var modelData
                    Layout.fillWidth: true
                    implicitHeight: 24
                    radius: 3
                    color: !!(root.win.activeCanvas && root.win.activeCanvas.currentArrowMode === amRow.modelData.mode) ? Theme.selected : (amRowMouse.containsMouse ? Theme.hover : "transparent")
                    Behavior on color { ColorAnimation { duration: 150 } }
                    Image {
                        id: amRowIcon
                        anchors { left: parent ? parent.left : undefined; leftMargin: 6; verticalCenter: parent ? parent.verticalCenter : undefined }
                        source: "icons/" + amRow.modelData.icon
                        width: 18
                        height: 18
                        sourceSize: Qt.size(18, 18)
                        onStatusChanged: {
                            if (status === Image.Error) {
                                console.warn("Failed to load icon:", amRow.modelData.icon)
                            }
                        }
                    }
                    Text {
                        anchors { left: amRowIcon.right; leftMargin: 6; verticalCenter: parent ? parent.verticalCenter : undefined }
                        text: amRow.modelData.label
                        color: !!(root.win.activeCanvas && root.win.activeCanvas.currentArrowMode === amRow.modelData.mode) ? Theme.accent : Theme.textPrimary
                        font { pixelSize: Theme.fontSizeLabel; family: Theme.fontDisplay }
                    }
                    MouseArea {
                        id: amRowMouse
                        anchors.fill: parent
                        hoverEnabled: true
                        cursorShape: Qt.PointingHandCursor
                        onClicked: {
                            if (root.win.activeCanvas) root.win.activeCanvas.currentArrowMode = amRow.modelData.mode
                            reactionsMenuTrigger.close()
                        }
                    }
                }
            }

            Rectangle { Layout.fillWidth: true; Layout.preferredHeight: 1; color: Theme.outline }
            Text {
                text: "Atom mapping"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontDisplay }
            }
            MenuItemRow {
                Layout.fillWidth: true
                text: "Auto-map Reaction"
                iconSource: "reaction-map.svg"
                onTriggered: {
                    root.win.isProcessing = true
                    if (root.win.activeSketch) root.win.activeSketch.requestSerialize("automap")
                    reactionsMenuTrigger.close()
                }
            }
            MenuItemRow {
                Layout.fillWidth: true
                text: "Clear Mapping"
                iconSource: "reaction-map.svg"
                onTriggered: {
                    root.win.isProcessing = true
                    if (root.win.activeSketch) root.win.activeSketch.requestSerialize("clear_mapping")
                    reactionsMenuTrigger.close()
                }
            }
            MenuItemRow {
                Layout.fillWidth: true
                text: "Correct Reacting Centers"
                iconSource: "reaction-map.svg"
                onTriggered: {
                    root.win.isProcessing = true
                    if (root.win.activeSketch) root.win.activeSketch.requestSerialize("correct_reacting_centers")
                    reactionsMenuTrigger.close()
                }
            }
        }
    }
    AppMenuBarItem {
        id: queryMenuTrigger
        text: "Query"
        ColumnLayout {
            width: 210
            spacing: 6

            Text {
                text: "Query atoms"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontDisplay }
            }
            GridLayout {
                columns: 4
                columnSpacing: Theme.spacingSmall
                rowSpacing: Theme.spacingSmall
                Layout.alignment: Qt.AlignHCenter
                Repeater {
                    model: [
                        { id: "ATOM_A",  icon: "A",  tip: "Any atom" },
                        { id: "ATOM_AH", icon: "AH", tip: "Any atom including H" },
                        { id: "ATOM_Q",  icon: "Q",  tip: "Any heteroatom" },
                        { id: "ATOM_QH", icon: "QH", tip: "Heteroatom or H" },
                        { id: "ATOM_M",  icon: "M",  tip: "Any metal" },
                        { id: "ATOM_MH", icon: "MH", tip: "Metal or H" },
                        { id: "ATOM_X",  icon: "X",  tip: "Any halogen" },
                        { id: "ATOM_XH", icon: "XH", tip: "Halogen or H" }
                    ]
                    delegate: IconCell {
                        id: qaCell
                        required property var modelData
                        cellWidth: 42
                        cellHeight: 32
                        glyph: qaCell.modelData.icon
                        tip: qaCell.modelData.tip
                        selected: !!(root.win.activeCanvas && root.win.activeCanvas.currentTool === qaCell.modelData.id)
                        onClicked: {
                            if (!root.win.activeCanvas) return
                            root.win.activeCanvas.currentTool = root.win.activeCanvas.currentTool === qaCell.modelData.id ? "SELECT" : qaCell.modelData.id
                            queryMenuTrigger.close()
                        }
                    }
                }
            }

            Rectangle { Layout.fillWidth: true; Layout.preferredHeight: 1; color: Theme.outline }

            Text {
                text: "R-group labels"
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeCaption; bold: true; letterSpacing: 1; family: Theme.fontDisplay }
            }
            GridLayout {
                columns: 4
                columnSpacing: Theme.spacingSmall
                rowSpacing: Theme.spacingSmall
                Layout.alignment: Qt.AlignHCenter
                Repeater {
                    model: ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]
                    delegate: IconCell {
                        id: rgCell
                        required property string modelData
                        cellWidth: 42
                        cellHeight: 28
                        glyph: rgCell.modelData
                        hasGlyphColorOverride: true
                        glyphColor: Theme.rgroupColor
                        accentColor: Theme.rgroupColor
                        borderAlwaysVisible: true
                        onClicked: {
                            if (root.win.activeCanvas) root.win.activeCanvas.currentTool = "ATOM_" + rgCell.modelData
                            queryMenuTrigger.close()
                        }
                    }
                }
            }
        }
    }
    AppMenuBarItem {
        text: "Copy"
        MenuItemRow { text: "Copy SMILES"; onTriggered: if (root.win.activeSketch) root.win.activeSketch.requestSerialize("smiles") }
        MenuItemRow { text: "Copy Canonical SMILES"; onTriggered: if (root.win.activeSketch) root.win.activeSketch.requestSerialize("canonical_smiles") }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Copy InChI"; onTriggered: if (root.win.activeSketch) root.win.activeSketch.requestSerialize("inchi") }
        MenuItemRow { text: "Copy InChIKey"; onTriggered: if (root.win.activeSketch) root.win.activeSketch.requestSerialize("inchikey") }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Copy Hash"; onTriggered: if (root.win.activeSketch) root.win.activeSketch.requestSerialize("hash") }
        MenuItemRow { text: "Copy Mass Composition"; onTriggered: if (root.win.activeSketch) root.win.activeSketch.requestSerialize("mass_composition") }
        MenuItemRow { text: "Copy pKa Values"; onTriggered: if (root.win.activeSketch) root.win.activeSketch.requestSerialize("pka_values") }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Compare Similarity…"; onTriggered: root.win.similarityDialog.open() }
    }
    AppMenuBarItem {
        text: "View"
        MenuItemRow { text: "Fit to Screen"; iconSource: "fit.svg"; shortcutHint: "Ctrl+0"; onTriggered: root.win.activeCanvas.fitToMolecule() }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { text: "Zoom In"; iconSource: "zoom-in.svg"; onTriggered: root.win.zoomLevel = Math.min(3.0, root.win.zoomLevel + 0.1) }
        MenuItemRow { text: "Zoom Out"; iconSource: "zoom-out.svg"; onTriggered: root.win.zoomLevel = Math.max(0.1, root.win.zoomLevel - 0.1) }
        MenuItemRow { text: "Reset (100%)"; onTriggered: root.win.zoomLevel = 1.0 }
        Rectangle { width: parent.width; height: 1; color: Theme.outline; opacity: 0.6 }
        MenuItemRow { 
            text: (root.win.uiSettings.cpkColorsEnabled ? "☑ " : "☐ ") + "CPK Atom Colors"
            onTriggered: {
                root.win.uiSettings.cpkColorsEnabled = !root.win.uiSettings.cpkColorsEnabled
                Theme.cpkColorsEnabled = root.win.uiSettings.cpkColorsEnabled
            }
        }
    }
    AppMenuBarItem {
        text: "Settings"
        MenuItemRow { text: "Drawing Settings…"; onTriggered: root.win.settingsDialogsGroup.settingsDialog.open() }
    }
    AppMenuBarItem {
        text: "Help"
        MenuItemRow { text: "Keyboard Shortcuts"; onTriggered: root.win.taskDialogsGroup.shortcutsDialog.open() }
    }
}
