import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import QtCore
import "../js/Selection.js" as Selection
import "../components"
import ".."

Item {
    id: root
    // Declared as a direct child of the root ApplicationWindow (real size), but
    // a plain Item never auto-fills its parent -- without this, root stayed at
    // its default 0x0 size at (0,0), so every FileDialog inside (anchors.centerIn:
    // parent, from components/TaskDialog.qml's shared base) centered on that
    // single point instead of the real window, rendering pinned to the
    // top-left corner with most of the dialog off-screen. Same bug class as
    // components/StatusPlaceholder.qml's earlier sizing fix.
    anchors.fill: parent
    required property var win

    property alias openDialog: openDialog
    property alias saveDialog: saveDialog
    property alias pdfSaveDialog: pdfSaveDialog
    property alias gridSaveDialog: gridSaveDialog
    property alias batchGridSaveDialog: batchGridSaveDialog
    property alias batchFileSaveDialog: batchFileSaveDialog
    property alias imageFileDialog: imageFileDialog

    FileDialog {
        id: openDialog
        title: "Open Molecule"
        fileMode: FileDialog.OpenFile
        nameFilters: [
            "Molfile (*.mol)", "SDF (*.sdf)", "RDF (*.rdf)", "SMILES (*.smi *.smiles)",
            "CML (*.cml)", "CDX (*.cdx)", "Ketcher JSON (*.ket)",
            "FASTA (*.fasta *.fa)", "HELM (*.helm)", "IDT Oligo (*.idt)",
            "All files (*)"
        ]
        onAccepted: win.loadFromFile(selectedFile)
    }

    FileDialog {
        id: saveDialog
        title: "Save Molecule"
        fileMode: FileDialog.SaveFile
        nameFilters: ["QtSketch Project (*.qtsketch)", "Molfile (*.mol)", "SDF (*.sdf)", "Ketcher JSON (*.ket)", "PNG image (*.png)", "SVG image (*.svg)", "All files (*)"]
        onAccepted: {
            const fileStr = selectedFile.toString().toLowerCase()
            if (fileStr.endsWith(".png")) {
                win.pendingRenderUrl = selectedFile
                // Full document snapshot (not a bare molfile) so text/images/
                // reaction arrows/etc. reach the exported file too.
                if (win.activeSketch) win.activeSketch.requestProjectSerialization("render_png")
                return
            }
            if (fileStr.endsWith(".svg")) {
                win.pendingRenderUrl = selectedFile
                if (win.activeSketch) win.activeSketch.requestProjectSerialization("render_svg")
                return
            }
            let fmt = "mol"
            if (fileStr.endsWith(".sdf")) fmt = "sdf"
            else if (fileStr.endsWith(".ket")) fmt = "ket"
            else if (fileStr.endsWith(".qtsketch")) fmt = "qtsketch"

            if (fmt === "qtsketch") {
                // Synchronous call with a real success/failure result -- only
                // record the filename and clear dirty state on confirmed
                // success. Never fall back to MOL on failure: that would
                // silently replace an existing native project with
                // chemistry-only content at the same destination.
                if (win.activeSketch) {
                    const success = win.activeSketch.saveProject(selectedFile)
                    if (success) {
                        win.pendingSaveUrl = selectedFile
                        win.setDocFile(DocumentManager.activeDocId, selectedFile)
                        if (win.activeCanvas) win.activeCanvas.setClean()
                    } else {
                        win.messageDialogsGroup.workerErrorDialog.errorText =
                            "Save failed: could not write the project file to " + selectedFile + ". The document was not modified on disk; your previous save (if any) is intact."
                        win.messageDialogsGroup.workerErrorDialog.severe = true
                        win.messageDialogsGroup.workerErrorDialog.open()
                    }
                }
            } else {
                win.pendingSaveUrl = selectedFile
                win.setDocFile(DocumentManager.activeDocId, selectedFile)
                if (win.activeSketch) win.activeSketch.requestStructure(fmt, "save")
            }
        }
    }

    FileDialog {
        id: pdfSaveDialog
        title: "Export as PDF"
        fileMode: FileDialog.SaveFile
        nameFilters: ["PDF document (*.pdf)"]
        onAccepted: {
            win.pendingRenderUrl = selectedFile
            if (win.activeSketch) win.activeSketch.requestProjectSerialization("render_pdf")
        }
    }

    FileDialog {
        id: gridSaveDialog
        title: "Export Reaction Scheme"
        fileMode: FileDialog.SaveFile
        nameFilters: ["PDF document (*.pdf)"]
        onAccepted: {
            win.pendingRenderUrl = selectedFile
            if (win.activeSketch) win.activeSketch.requestStructure("mol", "render_grid")
        }
    }

    FileDialog {
        id: batchGridSaveDialog
        title: "Export Batch as Image Grid"
        fileMode: FileDialog.SaveFile
        nameFilters: ["PDF document (*.pdf)", "PNG image (*.png)"]
        onAccepted: {
            win.pendingRenderUrl = selectedFile
            const fileStr = selectedFile.toString().toLowerCase()
            const format = fileStr.endsWith(".png") ? "png" : "pdf"
            win.indigoSvc.exportBatchGridToFile(win._pendingBatchGridMolfiles, selectedFile, format)
        }
    }

    FileDialog {
        id: batchFileSaveDialog
        title: "Export Batch to File"
        fileMode: FileDialog.SaveFile
        nameFilters: ["SDF (*.sdf)", "RDF (*.rdf)", "SMILES (*.smi)", "CML (*.cml)"]
        onAccepted: {
            const fileStr = selectedFile.toString().toLowerCase()
            let format = "sdf"
            if (fileStr.endsWith(".rdf")) format = "rdf"
            else if (fileStr.endsWith(".smi")) format = "smiles"
            else if (fileStr.endsWith(".cml")) format = "cml"
            win.indigoSvc.exportBatchToFile(win._pendingBatchGridMolfiles, selectedFile, format)
        }
    }

    FileDialog {
        id: imageFileDialog
        nameFilters: ["Images (*.png *.jpg *.jpeg *.gif *.bmp)"]

        property real chemX: 0
        property real chemY: 0
        property url pendingFileUrl
        // Tunable confidence threshold on Imago's reported recognition-warning count.
        // Verified empirically: a real chemical-structure image reported 0 warnings, an
        // unrelated screenshot reported 34 while still producing a (garbage) molfile — Imago
        // never refuses outright, so this count is the real signal to gate on. Not pinned at
        // exactly 0 since an imperfect scan may legitimately produce a few warnings and still
        // recognize correctly; revisit this number against real user images.
        property int maxAcceptableWarnings: 5

        onAccepted: {
            if (!win.activeSketch) return
            pendingFileUrl = selectedFile
            win.isProcessing = true
            win.imagoSvc.recognizeImage(selectedFile)
        }
    }

}
