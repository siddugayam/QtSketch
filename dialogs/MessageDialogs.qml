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
    // See dialogs/FileDialogs.qml's identical comment -- without this, root
    // stayed at its default 0x0 size at (0,0), pinning every MessageDialog
    // inside to the window's top-left corner instead of its real center.
    anchors.fill: parent
    required property var win

    property alias unsavedChangesDialog: unsavedChangesDialog
    property alias workerErrorDialog: workerErrorDialog
    property alias similarityResultDialog: similarityResultDialog
    property alias similarityRankResultDialog: similarityRankResultDialog

    MessageDialog {
        id: unsavedChangesDialog
        title: "Unsaved Changes"
        buttons: MessageDialog.Save | MessageDialog.Discard | MessageDialog.Cancel
        text: win._pendingCloseDocId >= 0 ? "This document has unsaved changes. Save before closing?" : "You have unsaved changes. Do you want to save them before exiting?"

        onButtonClicked: function(button) {
            if (button === MessageDialog.Save) {
                win.fileDialogsGroup.saveDialog.open()
            } else if (button === MessageDialog.Discard) {
                if (win._pendingCloseDocId >= 0) {
                    var docId = win._pendingCloseDocId
                    win._pendingCloseDocId = -1
                    // Look up the recovery id BEFORE closeDocument() removes it from
                    // DocumentManager's bookkeeping.
                    var rid = DocumentManager.recoveryIdFor(docId)
                    DocumentManager.closeDocument(docId)
                    const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
                    // Remove both new and old formats
                    win.fileIO.remove(appData + "/recovery/autosave_" + rid + ".qtsketch")
                    win.fileIO.remove(appData + "/recovery/autosave_" + rid + ".ket")
                } else {
                    if (win.discardRecoveryFiles) win.discardRecoveryFiles()
                    win.removeOpenAutosaveFiles()
                    win._forceQuit = true
                    Qt.quit()
                }
            }
            // Cancel: reset pending state
            if (button === MessageDialog.Cancel) {
                win._pendingCloseDocId = -1
            }
        }
    }

    MessageDialog {
        id: workerErrorDialog
        property string errorText: ""
        // Defaults to true (the scary "engine may have stopped" framing) since
        // that's correct for its real crash-signal caller (window's
        // onErrorOccurred, a genuine worker exception) -- MainWindow.qml's
        // ~14 other call sites for routine, recoverable feature failures
        // (batch/similarity/scaffold/render errors etc.) explicitly set this
        // to false right before opening. Deriving severity from errorText
        // keywords instead was tried and rejected: the real crash-path
        // messages ("Chemistry engine script not found...", worker JS
        // exception text) don't reliably contain any particular keyword, so
        // a keyword heuristic misclassified the actual severe case as mild.
        property bool severe: true
        title: severe ? "Chemistry Engine Error" : "Operation Failed"
        buttons: MessageDialog.Ok
        text: severe ?
              ("The chemistry engine encountered a severe error and may have stopped.\n\n" + errorText + "\n\nPlease save your work to a new file and restart the application.") :
              ("The operation could not be completed:\n\n" + errorText)
    }

    MessageDialog {
        id: similarityResultDialog
        title: "Similarity Result"
        text: ""
    }

    MessageDialog {
        id: similarityRankResultDialog
        title: "Similarity Ranking"
        text: ""
    }

    property alias recoveryDialog: recoveryDialog

    MessageDialog {
        id: recoveryDialog
        title: "Recover Unsaved Documents?"
        buttons: MessageDialog.Yes | MessageDialog.No
        text: "The application closed unexpectedly. Recover unsaved changes?"
        onButtonClicked: function(button) {
            if (button === MessageDialog.Yes) {
                win.restoreRecoveryFiles()
            } else if (button === MessageDialog.No) {
                win.discardRecoveryFiles()
            }
        }
    }

}
