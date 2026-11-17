import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
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
