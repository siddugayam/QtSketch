import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Popup {
    id: root
    width: 400
    height: Math.min(600, 60 + listView.contentHeight)
    x: parent ? (parent.width - width) / 2 : 0
    y: parent ? (parent.height - height) / 4 : 0
    modal: true
    focus: true
    closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutside
    
    required property var window

    property var commandRegistry: [
        { name: "File: New Document", action: function() { DocumentManager.addDocument() } },
        { name: "File: Open…", action: function() { window.openDialog.open() } },
        { name: "File: Save", action: function() { window.saveActive(false) } },
        { name: "File: Save As…", action: function() { window.saveActive(true) } },
        { name: "File: Load from SMILES…", action: function() { window.smilesDialog.open() } },
        { name: "File: Load from InChI…", action: function() { window.inchiLoadDialog.open() } },
        { name: "File: Biopolymer…", action: function() { window.biopolymerDialog.open() } },
        { name: "File: Export as PDF…", action: function() { window.printToPdf() } },
        { name: "File: Export Reaction Scheme (Grid)…", action: function() { window.gridSaveDialog.open() } },
        { name: "File: Clear Canvas", action: function() { if (window.activeCanvas) window.activeCanvas.clearCanvas() } },
        
        { name: "Edit: Undo", action: function() { if (window.activeCanvas && window.activeCanvas.canUndo) window.activeCanvas.undo() } },
        { name: "Edit: Redo", action: function() { if (window.activeCanvas && window.activeCanvas.canRedo) window.activeCanvas.redo() } },
        { name: "Edit: Cut", action: function() { if (window.activeCanvas) window.activeCanvas.cutSelection() } },
        { name: "Edit: Copy", action: function() { if (window.activeCanvas) window.activeCanvas.copySelection() } },
        { name: "Edit: Paste", action: function() { if (window.activeCanvas) window.activeCanvas.pasteSelection() } },
        { name: "Edit: Copy as Image", action: function() { if (window.activeCanvas) window.activeCanvas.copyAsImage() } },
        
        { name: "Structure: Layout", action: function() { window.executeStructureOp("layout") } },
        { name: "Structure: Layout Selected", action: function() { window.executeStructureOp("layoutSelected") } },
        { name: "Structure: Clean 2D", action: function() { window.executeStructureOp("clean2d") } },
        { name: "Structure: Aromatize", action: function() { window.executeStructureOp("aromatize") } },
        { name: "Structure: Dearomatize", action: function() { window.executeStructureOp("dearomatize") } },
        { name: "Structure: Add Explicit H", action: function() { window.executeStructureOp("unfoldH") } },
        { name: "Structure: Remove Explicit H", action: function() { window.executeStructureOp("foldH") } },
        { name: "Structure: Normalize", action: function() { window.executeStructureOp("normalize") } },
        { name: "Structure: Standardize", action: function() { window.executeStructureOp("standardize") } },
        { name: "Structure: Ionize at pH…", action: function() { window.ionizeDialog.open() } },
        { name: "Structure: Validate", action: function() { window.executeStructureOp("check") } },
        { name: "Structure: Search Substructure (SMARTS)…", action: function() { window.executeStructureOp("smartsSearch") } },
        { name: "Structure: R-Groups…", action: function() { window.executeStructureOp("rgroups") } },
        
        { name: "Reaction: Auto-map Reaction", action: function() { window.isProcessing = true; if (window.activeSketch) window.activeSketch.requestSerialize("automap") } },
        { name: "Reaction: Clear Mapping", action: function() { window.isProcessing = true; if (window.activeSketch) window.activeSketch.requestSerialize("clear_mapping") } },
        { name: "Reaction: Correct Reacting Centers", action: function() { window.isProcessing = true; if (window.activeSketch) window.activeSketch.requestSerialize("correct_reacting_centers") } },
        
        { name: "Copy: SMILES", action: function() { if (window.activeSketch) window.activeSketch.requestSerialize("smiles") } },
        { name: "Copy: Canonical SMILES", action: function() { if (window.activeSketch) window.activeSketch.requestSerialize("canonical_smiles") } },
        { name: "Copy: InChI", action: function() { if (window.activeSketch) window.activeSketch.requestSerialize("inchi") } },
        { name: "Copy: InChIKey", action: function() { if (window.activeSketch) window.activeSketch.requestSerialize("inchikey") } },
        { name: "Copy: Hash", action: function() { if (window.activeSketch) window.activeSketch.requestSerialize("hash") } },
        { name: "Copy: Mass Composition", action: function() { if (window.activeSketch) window.activeSketch.requestSerialize("mass_composition") } },
        { name: "Copy: pKa Values", action: function() { if (window.activeSketch) window.activeSketch.requestSerialize("pka_values") } },
        { name: "Copy: Compare Similarity…", action: function() { window.similarityDialog.open() } },
        
        { name: "View: Fit to Screen", action: function() { if (window.activeCanvas) window.activeCanvas.fitToMolecule() } },
        { name: "View: Zoom In", action: function() { window.zoomLevel = Math.min(3.0, window.zoomLevel + 0.1) } },
        { name: "View: Zoom Out", action: function() { window.zoomLevel = Math.max(0.1, window.zoomLevel - 0.1) } },
        { name: "View: Reset (100%)", action: function() { window.zoomLevel = 1.0 } },
        { name: "View: Toggle CPK Atom Colors", action: function() { 
            window.uiSettings.cpkColorsEnabled = !window.uiSettings.cpkColorsEnabled; 
            Theme.cpkColorsEnabled = window.uiSettings.cpkColorsEnabled; 
            if (window.activeCanvas) window.activeCanvas.refresh() 
        } }
    ]

    property var filteredCommands: commandRegistry

    onOpened: {
        searchInput.text = ""
        searchInput.forceActiveFocus()
        filteredCommands = commandRegistry
        listView.currentIndex = 0
    }

    function executeSelected() {
        if (listView.currentIndex >= 0 && listView.currentIndex < filteredCommands.length) {
            var cmd = filteredCommands[listView.currentIndex]
            root.close()
            cmd.action()
        }
    }

    background: Rectangle {
        color: Theme.surface
        border.color: Theme.outline
        radius: 8
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 10
        spacing: 10

        TextField {
            id: searchInput
            Layout.fillWidth: true
            placeholderText: "Search commands..."
            font { pixelSize: Theme.fontSizeBody; family: Theme.fontFamily }
            
            onTextEdited: {
                var query = text.toLowerCase()
                var res = []
                for (var i = 0; i < commandRegistry.length; i++) {
                    if (commandRegistry[i].name.toLowerCase().indexOf(query) !== -1) {
                        res.push(commandRegistry[i])
                    }
                }
                filteredCommands = res
                listView.currentIndex = res.length > 0 ? 0 : -1
            }

            Keys.onUpPressed: {
                if (listView.currentIndex > 0) listView.currentIndex--
            }
            Keys.onDownPressed: {
                if (listView.currentIndex < filteredCommands.length - 1) listView.currentIndex++
            }
            Keys.onReturnPressed: {
                root.executeSelected()
            }
            Keys.onEnterPressed: {
                root.executeSelected()
            }
        }

        ListView {
            id: listView
            Layout.fillWidth: true
            Layout.fillHeight: true
            clip: true
            model: root.filteredCommands
            
            delegate: Rectangle {
                width: ListView.view.width
                height: 36
                color: ListView.isCurrentItem ? Theme.selected : (mouseArea.containsMouse ? Theme.hover : "transparent")
                radius: 4

                Text {
                    anchors.verticalCenter: parent.verticalCenter
                    anchors.left: parent.left
                    anchors.leftMargin: 10
                    text: modelData.name
                    color: (ListView.isCurrentItem || mouseArea.containsMouse) ? Theme.accent : Theme.textPrimary
                    font { pixelSize: Theme.fontSizeBody; family: Theme.fontDisplay }
                }

                MouseArea {
                    id: mouseArea
                    anchors.fill: parent
                    hoverEnabled: true
                    onClicked: {
                        listView.currentIndex = index
                        root.executeSelected()
                    }
                }
            }
        }
    }
}
