import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

TaskDialog {
    id: root
    title: "Biopolymer Editor"
    width: 720
    standardButtons: Dialog.Close

    // Emitted by Load button — caller calls the relevant IndigoService method
    signal loadRequested(string format, string seqType, string text)
    // Emitted by "Preview as Sequence View" — caller builds the native monomer-box view
    // (plain Sequence/FASTA text only, not exploded to atoms)
    signal previewSequenceViewRequested(string seqType, string text)
    // Emitted by Export buttons — caller calls requestSerialize with the right reqId
    signal exportRequested(string format)

    // Called by MainWindow when an export result arrives
    function showExportResult(label, text) {
        if (text && text.trim().length > 0) {
            bioText.text = text
            bioText.readOnly = true
            statusLabel.text = label + " ready — edit above or Copy All."
            statusLabel.color = Theme.textSecondary
        } else {
            statusLabel.text = "Export failed: structure may not be a biopolymer."
            statusLabel.color = "#c0392b"
        }
    }

    function showLoadError(error) {
        statusLabel.text = "Load error: " + error
        statusLabel.color = "#c0392b"
    }

    function onLoadSuccess() {
        statusLabel.text = "Loaded successfully."
        statusLabel.color = "#27ae60"
        root.close()
    }

    // Format / type configuration
    readonly property var formats:  ["Sequence", "FASTA", "HELM", "IDT", "AxoLabs", "SVG"]
    readonly property var seqTypes: ["DNA", "RNA", "Peptide"]
    
    // Vector export formats
    readonly property var vectorFormats: ["SVG", "PDF"]

    // IDT and AxoLabs are oligonucleotide-only — disable Peptide type for them
    property bool typeRelevant: formatCombo.currentIndex <= 1  // Sequence or FASTA
    property bool peptideDisabledForFormat: formatCombo.currentIndex >= 3  // IDT/AxoLabs

    // Reset textarea whenever format changes so stale data is not accidentally loaded
    onVisibleChanged: if (visible) _resetInput()

    function _resetInput() {
        bioText.text = ""
        bioText.readOnly = false
        statusLabel.text = ""
        // Auto-fix type when IDT/AxoLabs selected and Peptide was chosen
        if (peptideDisabledForFormat && seqTypeCombo.currentIndex === 2)
            seqTypeCombo.currentIndex = 0
    }

    ColumnLayout {
        width: root.availableWidth
        spacing: 8

        // ── Format & type row ───────────────────────────────────────────────
        RowLayout {
            spacing: 12
            Layout.fillWidth: true

            Label { text: "Format:"; color: Theme.textSecondary; font.pixelSize: 12 }
            ComboBox {
                id: formatCombo
                model: root.formats
                Layout.preferredWidth: 120
                onCurrentIndexChanged: root._resetInput()
            }

            Label {
                text: "Type:"
                color: root.typeRelevant ? Theme.textSecondary : Theme.textSecondary
                opacity: root.typeRelevant ? 1.0 : 0.35
                font.pixelSize: 12
            }
            ComboBox {
                id: seqTypeCombo
                model: root.seqTypes
                Layout.preferredWidth: 110
                enabled: root.typeRelevant
                opacity: root.typeRelevant ? 1.0 : 0.35
                // Remove Peptide option for IDT/AxoLabs
                onEnabledChanged: {
                    if (!enabled && currentIndex === 2) currentIndex = 0
                }
            }
        }

        // ── Text area ───────────────────────────────────────────────────────
        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 220
            color: bioText.readOnly ? Qt.darker(Theme.background, 1.05) : Theme.background
            border.color: Theme.outline
            radius: 4
            clip: true

            ScrollView {
                anchors { fill: parent; margins: 2 }
                contentWidth: availableWidth

                TextArea {
                    id: bioText
                    width: parent.width
                    font { family: Theme.fontMono; pixelSize: 12 }
                    wrapMode: TextArea.Wrap
                    selectByMouse: true
                    background: null
                    color: Theme.textPrimary
                    placeholderText: {
                        const fmt = root.formats[formatCombo.currentIndex]
                        const typ = root.seqTypes[seqTypeCombo.currentIndex]
                        if (fmt === "Sequence") {
                            if (typ === "DNA")     return "e.g.  ATGCATGCTAGC"
                            if (typ === "RNA")     return "e.g.  AUGCAUGCUAGC"
                            return "e.g.  ACDEFGHIKLMNPQRSTVWY"
                        }
                        if (fmt === "FASTA") return ">MySeq\nATGCATGCTAGC"
                        if (fmt === "HELM")  return "RNA1{R(A)P.R(U)P.R(G)P.R(C)}$$$$"
                        if (fmt === "IDT")   return "e.g.  /5Biosg/ATGCATGC/3AmMO/"
                        if (fmt === "AxoLabs") return "AxoLabs oligo notation"
                        return ""
                    }
                }
            }
        }

        // ── Status label ─────────────────────────────────────────────────────
        Label {
            id: statusLabel
            text: ""
            font.pixelSize: 11
            Layout.fillWidth: true
            wrapMode: Text.WordWrap
        }

        // ── Action buttons ────────────────────────────────────────────────────
        RowLayout {
            spacing: 6
            Layout.fillWidth: true

            Button {
                text: "Load to Canvas"
                enabled: bioText.text.trim().length > 0 && !bioText.readOnly
                onClicked: {
                    const fmt     = root.formats[formatCombo.currentIndex]
                    const seqType = root.seqTypes[seqTypeCombo.currentIndex].toUpperCase()
                    // Set "Loading…" BEFORE emitting: some validation failures (e.g. an
                    // unsafe sequence alphabet) are reported synchronously via
                    // biopolymerLoadError before loadRequested() even returns, so
                    // showLoadError() must run first and this line must not clobber it.
                    statusLabel.text  = "Loading…"
                    statusLabel.color = Theme.textSecondary
                    root.loadRequested(fmt, seqType, bioText.text.trim())
                }
            }

            Button {
                text: "Preview as Sequence View"
                visible: root.typeRelevant
                enabled: bioText.text.trim().length > 0 && !bioText.readOnly
                onClicked: {
                    const seqType = root.seqTypes[seqTypeCombo.currentIndex].toUpperCase()
                    root.previewSequenceViewRequested(seqType, bioText.text.trim())
                }
            }

            Button {
                text: "Clear"
                onClicked: root._resetInput()
            }

            Item { Layout.fillWidth: true }

            Button {
                text: "Copy All"
                enabled: bioText.text.trim().length > 0
                onClicked: {
                    bioText.selectAll()
                    bioText.copy()
                    bioText.deselect()
                    statusLabel.text = "Copied to clipboard."
                    statusLabel.color = Theme.textSecondary
                }
            }
        }

        RowLayout {
            spacing: 6
            Layout.fillWidth: true

            Label {
                text: "Export:"
                color: Theme.textSecondary
                font.pixelSize: 12
            }

            // Export buttons
            Repeater {
                model: root.formats
                Button {
                    required property string modelData
                    required property int index
                    text: modelData
                    implicitWidth: Math.max(68, implicitContentWidth + leftPadding + rightPadding)
                    // SVG is always enabled for vector export
                    opacity: (modelData === "SVG" || modelData === "PDF") ? 1.0 : 
                            (index >= 3 && seqTypeCombo.currentIndex === 2) ? 0.35 : 1.0
                    enabled: opacity > 0.5
                    onClicked: {
                        if (modelData === "SVG" || modelData === "PDF") {
                            // Vector export - requires document to be loaded
                            statusLabel.text = "Vector export requires an open document"
                            statusLabel.color = Theme.caution
                        } else {
                            bioText.readOnly = false
                            bioText.text = ""
                            statusLabel.text = "Exporting as " + modelData + "…"
                            statusLabel.color = Theme.textSecondary
                            root.exportRequested(modelData)
                        }
                    }
                }
            }

            Item { Layout.fillWidth: true }
        }
    }
}
