import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import Sketch.App

// Read-only "snake view" of a biopolymer as connected monomer boxes, instead of the
// exploded-atom structure the existing Load-to-Canvas path produces. Fully independent
// of chem-core.js's own (non-functional in this bundle -- see Features To Be Implemented.md /
// Architecture Map.md, 2026-07-23 entries) DrawingEntitiesManager engine: the worker just
// hands back plain {id,label,x,y,monomerClass} monomers and {fromId,toId} bonds, painted here.
TaskDialog {
    id: root
    title: "Sequence View"
    width: 640
    height: 420
    standardButtons: Dialog.Close

    property var snapshot: null
    property int selectedMonomerId: -1
    property var sketch: null

    FileIO { id: svgFileIO }

    FileDialog {
        id: svgSaveDialog
        title: "Save SVG"
        fileMode: FileDialog.SaveFile
        nameFilters: ["SVG image (*.svg)", "All files (*)"]
        onAccepted: {
            const svg = root.exportAsSvg()
            if (svg && svgFileIO.write(selectedFile, svg)) {
                statusLabel.text = "SVG saved to " + selectedFile
                statusLabel.color = Theme.success
            } else {
                statusLabel.text = "Save failed"
                statusLabel.color = Theme.error
            }
        }
    }

    function openSnapshot(dataJson) {
        try {
            root.snapshot = JSON.parse(dataJson)
        } catch (e) {
            root.snapshot = null
        }
        root.selectedMonomerId = -1
        root.open()
        seqCanvas.requestPaint()
    }

    onVisibleChanged: if (visible) seqCanvas.requestPaint()

    readonly property int cellWidth: 60
    readonly property int cellHeight: 120
    readonly property int boxSize: 36

    function monomerColor(monomerClass) {
        return monomerClass === "Base" ? Theme.badgeReactingCenter : Theme.accent
    }

    function monomerById(id) {
        if (!root.snapshot || !root.snapshot.monomers) return null
        for (let i = 0; i < root.snapshot.monomers.length; i++) {
            if (root.snapshot.monomers[i].id === id) return root.snapshot.monomers[i]
        }
        return null
    }

    // Export the current view as SVG
    function exportAsSvg() {
        if (!root.snapshot || !root.snapshot.monomers) return ""
        
        const originX = 24, originY = 24
        const half = root.boxSize / 2
        
        // Calculate bounds
        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity
        for (let i = 0; i < root.snapshot.monomers.length; i++) {
            const m = root.snapshot.monomers[i]
            minX = Math.min(minX, m.x)
            maxX = Math.max(maxX, m.x + root.boxSize)
            minY = Math.min(minY, m.y)
            maxY = Math.max(maxY, m.y + root.boxSize)
        }
        
        const width = Math.max(800, (maxX - minX + 40) * 2)
        const height = Math.max(600, (maxY - minY + 40) * 2)
        
        let svg = '<?xml version="1.0" encoding="UTF-8"?>\n'
        svg += `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">\n`
        svg += '<defs/>\n'
        
        // Add connections
        const bonds = root.snapshot.bonds || []
        for (let i = 0; i < bonds.length; i++) {
            const a = root.monomerById(bonds[i].fromId)
            const b = root.monomerById(bonds[i].toId)
            if (!a || !b) continue
            const x1 = originX + a.x + half
            const y1 = originY + a.y + half
            const x2 = originX + b.x + half
            const y2 = originY + b.y + half
            svg += `<line x1="${x1 * 2}" y1="${y1 * 2}" x2="${x2 * 2}" y2="${y2 * 2}" stroke="#212121" stroke-width="2"/>`
        }
        
        // Add residues
        const monomers = root.snapshot.monomers
        for (let i = 0; i < monomers.length; i++) {
            const m = monomers[i]
            const x = originX + m.x
            const y = originY + m.y
            const color = root.monomerColor(m.monomerClass)
            const isSelected = m.id === root.selectedMonomerId
            
            // Box
            svg += `<rect x="${x * 2}" y="${y * 2}" width="${root.boxSize * 2}" height="${root.boxSize * 2}" fill="${color}" stroke="${isSelected ? '#000000' : Qt.darker(color, 1.3)}" stroke-width="${isSelected ? 3 : 1}" rx="6" ry="6"/>`
            
            // Label
            svg += `<text x="${(x + half / 2) * 2}" y="${(y + half) * 2}" font-family="Segoe UI" font-size="12" text-anchor="middle" dominant-baseline="central" fill="#FFFFFF">${m.label}</text>`
            
            // Number
            svg += `<text x="${(x + 2) * 2}" y="${(y + 10) * 2}" font-family="Segoe UI" font-size="10" text-anchor="start" dominant-baseline="hanging" fill="#000000">${i + 1}</text>`
        }
        
        // Add direction indicator
        if (monomers.length > 0) {
            const first = monomers[0]
            const last = monomers[monomers.length - 1]
            const isProtein = root.snapshot.seqType === 'PEPTIDE'
            const direction = isProtein ? 'N→C' : "5'→3'"
            
            svg += `<line x1="${(originX + first.x - 20) * 2}" y1="${(originY + first.y + half) * 2}" x2="${(originX + first.x - 5) * 2}" y2="${(originY + first.y + half) * 2}" stroke="#0066CC" stroke-width="2"/>`
            svg += `<text x="${(originX + first.x - 25) * 2}" y="${(originY + first.y + half - 5) * 2}" font-family="Segoe UI" font-size="10" text-anchor="end" dominant-baseline="central" fill="#0066CC">${direction}</text>`
        }
        
        svg += '</svg>'
        return svg
    }

    // Copy SVG to clipboard
    function copySvgToClipboard() {
        const svg = root.exportAsSvg()
        if (svg && root.sketch) {
            root.sketch.setOsClipboardText(svg)
            statusLabel.text = "SVG copied to clipboard"
            statusLabel.color = Theme.success
        }
    }

    // Save SVG to a file, via a real native save dialog
    function downloadSvgFile() {
        if (!root.exportAsSvg()) return
        svgSaveDialog.open()
    }

    ColumnLayout {
        anchors { fill: parent; margins: 8 }
        spacing: 6

        Label {
            text: root.snapshot && root.snapshot.monomers
                  ? root.snapshot.monomers.length + " monomer(s), " + root.snapshot.bonds.length + " bond(s) — click a monomer to inspect"
                  : "No sequence loaded."
            color: Theme.textSecondary
            font.pixelSize: 11
        }

        RowLayout {
            Layout.fillWidth: true
            spacing: 6

            TextField {
                id: addSymbolField
                Layout.preferredWidth: 60
                placeholderText: "e.g. X"
                maximumLength: 1
            }
            Button {
                text: "Add to End"
                enabled: root.sketch && addSymbolField.text.trim().length > 0
                onClicked: {
                    const seqType = (root.snapshot && root.snapshot.seqType) ? root.snapshot.seqType : "PEPTIDE"
                    root.sketch.sendCommand("bioAddMonomer", [addSymbolField.text.trim(), seqType])
                    addSymbolField.text = ""
                }
            }
            Button {
                text: "Delete Selected"
                enabled: root.sketch && root.selectedMonomerId >= 0
                onClicked: {
                    root.sketch.sendCommand("bioDeleteMonomer", [root.selectedMonomerId])
                    root.selectedMonomerId = -1
                }
            }
            Item { Layout.fillWidth: true }

            // Vector export buttons
            Button {
                text: "Copy SVG"
                enabled: root.snapshot && root.snapshot.monomers && root.snapshot.monomers.length > 0
                onClicked: root.copySvgToClipboard()
                ToolTip.text: "Copy sequence view as SVG to clipboard"
                ToolTip.visible: hovered
            }
            Button {
                text: "Save SVG"
                enabled: root.snapshot && root.snapshot.monomers && root.snapshot.monomers.length > 0
                onClicked: root.downloadSvgFile()
                ToolTip.text: "Save sequence view as SVG file"
                ToolTip.visible: hovered
            }
        }

        // Status label for export operations
        Label {
            id: statusLabel
            text: ""
            color: Theme.textSecondary
            font.pixelSize: 11
            Layout.fillWidth: true
            wrapMode: Text.WordWrap
            visible: text.length > 0
        }

        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: Theme.background
            border.color: Theme.outline
            radius: 4
            clip: true

            Flickable {
                id: flick
                anchors.fill: parent
                anchors.margins: 1
                contentWidth: Math.max(width, seqCanvas.width)
                contentHeight: Math.max(height, seqCanvas.height)
                clip: true
                ScrollBar.horizontal: ScrollBar {}
                ScrollBar.vertical: ScrollBar {}

                Canvas {
                    id: seqCanvas
                    width: {
                        if (!root.snapshot || !root.snapshot.monomers || root.snapshot.monomers.length === 0) return flick.width
                        let maxX = 0
                        for (let i = 0; i < root.snapshot.monomers.length; i++)
                            maxX = Math.max(maxX, root.snapshot.monomers[i].x)
                        return Math.max(flick.width, maxX + root.cellWidth + 40)
                    }
                    height: {
                        if (!root.snapshot || !root.snapshot.monomers || root.snapshot.monomers.length === 0) return flick.height
                        let maxY = 0
                        for (let i = 0; i < root.snapshot.monomers.length; i++)
                            maxY = Math.max(maxY, root.snapshot.monomers[i].y)
                        return Math.max(flick.height, maxY + root.cellHeight)
                    }

                    onPaint: {
                        const ctx = getContext("2d")
                        ctx.reset()
                        if (!root.snapshot || !root.snapshot.monomers) return

                        const originX = 24, originY = 24
                        const half = root.boxSize / 2

                        // Bonds first, so boxes paint on top of the connecting lines
                        ctx.strokeStyle = Theme.outline
                        ctx.lineWidth = 2
                        const bonds = root.snapshot.bonds || []
                        for (let i = 0; i < bonds.length; i++) {
                            const a = root.monomerById(bonds[i].fromId)
                            const b = root.monomerById(bonds[i].toId)
                            if (!a || !b) continue
                            ctx.beginPath()
                            ctx.moveTo(originX + a.x + half, originY + a.y + half)
                            ctx.lineTo(originX + b.x + half, originY + b.y + half)
                            ctx.stroke()
                        }

                        const monomers = root.snapshot.monomers
                        for (let i = 0; i < monomers.length; i++) {
                            const m = monomers[i]
                            const x = originX + m.x, y = originY + m.y
                            ctx.fillStyle = root.monomerColor(m.monomerClass)
                            ctx.strokeStyle = (m.id === root.selectedMonomerId) ? Theme.textPrimary : Qt.darker(root.monomerColor(m.monomerClass), 1.3)
                            ctx.lineWidth = (m.id === root.selectedMonomerId) ? 3 : 1
                            ctx.beginPath()
                            if (ctx.roundRect) ctx.roundRect(x, y, root.boxSize, root.boxSize, 6)
                            else ctx.rect(x, y, root.boxSize, root.boxSize)
                            ctx.fill()
                            if (m.ambiguous) ctx.setLineDash([3, 2])
                            ctx.stroke()
                            if (m.ambiguous) ctx.setLineDash([])

                            ctx.fillStyle = Theme.badgeText
                            ctx.font = "bold 12px " + Theme.fontDisplayCss
                            ctx.textAlign = "center"
                            ctx.textBaseline = "middle"
                            ctx.fillText(m.label, x + half, y + half)
                        }
                    }

                    MouseArea {
                        anchors.fill: parent
                        onClicked: (m) => {
                            if (!root.snapshot || !root.snapshot.monomers) return
                            const originX = 24, originY = 24
                            let hit = -1
                            const monomers = root.snapshot.monomers
                            for (let i = 0; i < monomers.length; i++) {
                                const mm = monomers[i]
                                const x = originX + mm.x, y = originY + mm.y
                                if (m.x >= x && m.x <= x + root.boxSize && m.y >= y && m.y <= y + root.boxSize) {
                                    hit = mm.id
                                    break
                                }
                            }
                            root.selectedMonomerId = hit
                            seqCanvas.requestPaint()
                        }
                    }
                }
            }
        }

        Label {
            visible: root.selectedMonomerId >= 0
            text: {
                const m = root.monomerById(root.selectedMonomerId)
                return m ? ("Selected: " + m.alias + " (" + m.monomerClass + ")") : ""
            }
            color: Theme.textPrimary
            font.pixelSize: 11
        }
    }
}
