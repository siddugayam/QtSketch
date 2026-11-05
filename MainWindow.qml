import QtQuick
import QtQuick.Controls
import QtQuick.Templates as T
import QtQuick.Layouts
import QtQuick.Dialogs
import QtCore
import "js/Selection.js" as Selection
import "js/ToolLabels.js" as ToolLabels

ApplicationWindow {
    id: window
    property alias fileDialogsGroup: fileDialogsGroup
    property alias taskDialogsGroup: taskDialogsGroup
    property alias messageDialogsGroup: messageDialogsGroup

    flags: Qt.Window | Qt.FramelessWindowHint

    // Fluent's SmokeFillColorDefault: the scrim behind a modal content dialog.
    Overlay.modal: Rectangle { color: Theme.smokeFill }

    // Caption buttons to the Windows spec: 46x32 each, flush against one
    // another and the window's right edge, 10px Segoe Fluent Icons glyph. The
    // previous version used stock Buttons, which carry a 100x40 background in
    // every style — that width, not any spacing, is what pushed them apart.
    component CaptionButton: T.Button {
        id: capBtn
        property string glyph: ""
        property bool danger: false

        implicitWidth: Theme.captionButtonWidth
        implicitHeight: Theme.captionBarHeight
        focusPolicy: Qt.NoFocus

        background: Rectangle {
            color: capBtn.danger
                   ? (capBtn.down ? Qt.darker(Theme.closeButtonHover, 1.2)
                                  : (capBtn.hovered ? Theme.closeButtonHover : "transparent"))
                   : (capBtn.down ? Theme.pressed
                                  : (capBtn.hovered ? Theme.hover : "transparent"))
        }
        contentItem: Text {
            text: capBtn.glyph
            // Segoe Fluent Icons renders caption glyphs at 10px per the Windows
            // spec; they are drawn to look correct at that size, not scaled down.
            font { family: Theme.fontIcons; pixelSize: 10 }
            color: (capBtn.danger && (capBtn.hovered || capBtn.down)) ? "#FFFFFF" : Theme.textPrimary
            horizontalAlignment: Text.AlignHCenter
            verticalAlignment: Text.AlignVCenter
        }
    }

    header: Rectangle {
        height: Theme.captionBarHeight
        color: Theme.surface
        RowLayout {
            anchors.fill: parent
            spacing: 0
            Text {
                text: window.title
                color: Theme.textPrimary
                font { pixelSize: Theme.fontSizeCaption; family: Theme.fontFamily }
                Layout.leftMargin: 16
                Layout.fillWidth: true
                elide: Text.ElideRight
                verticalAlignment: Text.AlignVCenter
            }
            CaptionButton {
                glyph: ""                       // ChromeMinimize
                Accessible.name: "Minimize"
                onClicked: window.showMinimized()
            }
            CaptionButton {
                // ChromeRestore when maximized, ChromeMaximize otherwise
                glyph: window.visibility === Window.Maximized ? "" : ""
                Accessible.name: window.visibility === Window.Maximized ? "Restore down" : "Maximize"
                onClicked: {
                    if (window.visibility === Window.Maximized)
                        window.showNormal()
                    else
                        window.showMaximized()
                }
            }
            CaptionButton {
                glyph: ""                       // ChromeClose
                danger: true
                Accessible.name: "Close"
                onClicked: window.close()
            }
        }
        TapHandler {
            acceptedButtons: Qt.LeftButton
            onTapped: if (tapCount === 2) {
                if (window.visibility === Window.Maximized)
                    window.showNormal()
                else
                    window.showMaximized()
            }
        }
        DragHandler {
            onActiveChanged: if(active) window.startSystemMove()
        }
    }

    // Resize handlers for frameless window. DragHandler is not an Item — it has no
    // x/y/width/height/anchors of its own (Qt 6.11 QML type reference, verified) — so
    // each edge/corner hit-strip needs its own wrapping Item; the handler then acts
    // within that Item's bounds (its documented "parent" scope), not the handler's.
    Item {
        anchors.fill: parent
        Item {
            anchors { top: parent.top; left: parent.left; right: parent.right }
            height: 4
            DragHandler { target: null; cursorShape: Qt.SizeVerCursor; onActiveChanged: if (active) window.startSystemResize(Qt.TopEdge) }
        }
        Item {
            anchors { bottom: parent.bottom; left: parent.left; right: parent.right }
            height: 4
            DragHandler { target: null; cursorShape: Qt.SizeVerCursor; onActiveChanged: if (active) window.startSystemResize(Qt.BottomEdge) }
        }
        Item {
            anchors { left: parent.left; top: parent.top; bottom: parent.bottom }
            width: 4
            DragHandler { target: null; cursorShape: Qt.SizeHorCursor; onActiveChanged: if (active) window.startSystemResize(Qt.LeftEdge) }
        }
        Item {
            anchors { right: parent.right; top: parent.top; bottom: parent.bottom }
            width: 4
            DragHandler { target: null; cursorShape: Qt.SizeHorCursor; onActiveChanged: if (active) window.startSystemResize(Qt.RightEdge) }
        }
        Item {
            anchors { top: parent.top; left: parent.left }
            width: 6; height: 6
            DragHandler { target: null; cursorShape: Qt.SizeFDiagCursor; onActiveChanged: if (active) window.startSystemResize(Qt.TopLeftCorner) }
        }
        Item {
            anchors { top: parent.top; right: parent.right }
            width: 6; height: 6
            DragHandler { target: null; cursorShape: Qt.SizeBDiagCursor; onActiveChanged: if (active) window.startSystemResize(Qt.TopRightCorner) }
        }
        Item {
            anchors { bottom: parent.bottom; left: parent.left }
            width: 6; height: 6
            DragHandler { target: null; cursorShape: Qt.SizeBDiagCursor; onActiveChanged: if (active) window.startSystemResize(Qt.BottomLeftCorner) }
        }
        Item {
            anchors { bottom: parent.bottom; right: parent.right }
            width: 6; height: 6
            DragHandler { target: null; cursorShape: Qt.SizeFDiagCursor; onActiveChanged: if (active) window.startSystemResize(Qt.BottomRightCorner) }
        }
    }

    property bool isProcessing: false
    property bool explicitCheckPending: false
    property real canvasScaleX: 1.0
    property real canvasScaleY: 1.0
    property real zoomLevel: 1.0
    property real pixelsPerCm: 37.8
    property var pageMargins: ({ left: 2.0, right: 2.0, top: 2.0, bottom: 2.0 })

    // The currently-visible ChemCanvas instance and its bound document. Resolved once the
    // per-document Repeater below has instantiated at least one activeCanvas.
    property var activeCanvas: null
    readonly property var activeSketch: DocumentManager.documentFor(DocumentManager.activeDocId)
    property var dirtyDocs: ({})
    property int _pendingCloseDocId: -1

    onActiveCanvasChanged: {
        if (activeCanvas && activeCanvas.sketch) {
            AppController.setV8Process(activeCanvas.sketch)
            // Properties/stereo/check badges only refresh for the active document —
            // kick the debounce so a newly focused tab isn't stale until its next edit.
            propUpdateTimer.restart()
            activeCanvas.sketch.sendCommand("getClipboardPreview", [])
        }
    }

    // Per-document file paths + display titles (UI-only state). titleRev bumps
    // force TabButton text bindings to re-evaluate after open/save-as/rename.
    property var docFilePaths: ({})
    property int titleRev: 0
    property url pendingSaveUrl
    property url pendingRenderUrl
    function setDocFile(docId, fileUrl) {
        docFilePaths[docId] = fileUrl
        const s = fileUrl.toString()
        docTitles[docId] = decodeURIComponent(s.substring(s.lastIndexOf("/") + 1))
        titleRev++
        addRecentFile(fileUrl)
    }

    function addRecentFile(fileUrl) {
        if (!fileUrl) return
        var path = fileUrl.toString()
        if (path === "") return
        var files = uiSettings.recentFilesJoined ? uiSettings.recentFilesJoined.split("|") : []
        var idx = files.indexOf(path)
        if (idx !== -1) files.splice(idx, 1)
        files.unshift(path)
        if (files.length > 5) files = files.slice(0, 5)
        uiSettings.recentFilesJoined = files.join("|")
    }
    // Declared on window (root), not a nested Item, and called everywhere as
    // window.executeStructureOp(...) — the Structure menu's MenuItemRow children get
    // reparented into AppMenuBarItem's internal Popup (default property alias), and the
    // Qt Quick Compiler (qmlsc AOT) cannot resolve a bare, unqualified function name
    // declared on an un-id'd ancestor Item across that reparenting boundary (confirmed:
    // every OTHER menu in this same reparenting setup already calls its handlers via an
    // id-qualified path - window.saveActive(...)/window.printToPdf() in the File menu,
    // activeCanvas.xxx() in Edit - and those all work; this was the one unqualified
    // exception, and every Structure-menu item threw "ReferenceError: executeStructureOp
    // is not defined" at runtime as a result).
    function executeStructureOp(op) {
        if (op === "smartsSearch") { smartsSearchPopup.open(); return }
        if (op === "rgroups") { rgroupPanel.open(); return }
        if (!activeSketch) return
        if (op === "check") {
            window.explicitCheckPending = true
            activeSketch.requestSerialize("check")
            return
        }
        if (op === "layoutSelected") {
            activeSketch.sendCommand("layoutSelectedChain", [])
            return
        }
        if (op === "bracketSelection") {
            activeSketch.sendCommand("addBracketSelection", [])
            return
        }
        if (op === "iupac_name") {
            // Use special molfile that prevents fabricated E/Z labels
            activeSketch.requestMolfileForNaming("iupac_name")
            return
        }
        if (op !== "layout" && op !== "aromatize") window.isProcessing = true
        activeSketch.requestSerialize(op)
    }

    function saveActive(forceDialog) {
        if (!activeSketch) return
        const path = docFilePaths[DocumentManager.activeDocId]
        if (!forceDialog && path) {
            const fileStr = path.toString().toLowerCase()
            let fmt = "mol"
            if (fileStr.endsWith(".sdf")) fmt = "sdf"
            else if (fileStr.endsWith(".ket")) fmt = "ket"
            else if (fileStr.endsWith(".qtsketch")) fmt = "qtsketch"
            pendingSaveUrl = path
            addRecentFile(path)
            if (fmt === "qtsketch") {
                // Use new project save method for .qtsketch files
                const success = activeSketch.saveProject(path)
                if (!success) {
                    activeSketch.requestStructure("mol", "save")
                }
            } else {
                activeSketch.requestStructure(fmt, "save")
            }
            if (activeCanvas) activeCanvas.setClean()
            // Clean up recovery file on explicit save
            const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
            // Try to remove both new and old formats
            fileIO.remove(appData + "/recovery/autosave_" + DocumentManager.activeDocId + ".qtsketch")
            fileIO.remove(appData + "/recovery/autosave_" + DocumentManager.activeDocId + ".ket")
        } else {
            fileDialogsGroup.saveDialog.open()
        }
    }

    // Dialog ids below are declared deep in this file's own scope, which is
    // NOT automatically visible outside it (QML id-scoping is document-local —
    // it does not turn an id into a property of the containing root object).
    // AppMenus.qml holds a `win` reference to this window and reaches these
    // dialogs via e.g. `win.openDialog.open()`, so each one needed here MUST
    // be re-exposed as an alias, or that reference silently fails at runtime
    // (same failure class as the earlier executeStructureOp scoping bug).
    property alias openDialog: fileDialogsGroup.openDialog
    property alias gridSaveDialog: fileDialogsGroup.gridSaveDialog
    property alias biopolymerDialog: biopolymerDialog
    property alias smilesDialog: taskDialogsGroup.smilesDialog
    property alias similarityDialog: taskDialogsGroup.similarityDialog
    property alias ionizeDialog: taskDialogsGroup.ionizeDialog
    property alias inchiLoadDialog: taskDialogsGroup.inchiLoadDialog
    property alias rgroupPerMoleculeResultsDialog: taskDialogsGroup.rgroupPerMoleculeResultsDialog
    // Same reasoning: dialogs/FileDialogs.qml's batch export handlers reach
    // these two backend services via win.indigoSvc/win.imagoSvc — both need
    // the same alias treatment as the dialog ids above.
    property alias indigoSvc: indigoSvc
    property alias imagoSvc: imagoSvc
    property alias uiSettings: uiSettings

    property var docTitles: ({})

    function titleFor(docId) {
        if (docTitles[docId] === undefined) {
            var n = Object.keys(docTitles).length + 1
            docTitles[docId] = "Untitled " + n
        }
        return docTitles[docId]
    }

    property var _recoveredFiles: []

    function restoreRecoveryFiles() {
        const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
        // Sort files so .qtsketch files are restored first (newer format)
        const sortedFiles = window._recoveredFiles.sort((a, b) => {
            const aIsKet = a.endsWith(".ket")
            const bIsKet = b.endsWith(".ket")
            if (aIsKet && !bIsKet) return 1  // ket files last
            if (!aIsKet && bIsKet) return -1 // qtsketch files first
            return 0
        })
        for (let i = 0; i < sortedFiles.length; i++) {
            let path = appData + "/recovery/" + sortedFiles[i]
            window.loadFromFile(Qt.url("file:///" + path.replace(/\\/g, "/")))
        }
    }

    function discardRecoveryFiles() {
        const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
        for (let i = 0; i < window._recoveredFiles.length; i++) {
            fileIO.remove(appData + "/recovery/" + window._recoveredFiles[i])
        }
        window._recoveredFiles = []
    }

    property var _userTemplates: []

    function loadUserTemplates() {
        const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
        const path = appData + "/user_templates.json"
        if (fileIO.exists(path)) {
            const data = fileIO.read(path)
            try {
                window._userTemplates = JSON.parse(data)
            } catch (e) {
                window._userTemplates = []
            }
        }
    }
    
    function saveUserTemplates() {
        const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
        fileIO.mkpath(appData)
        fileIO.write(appData + "/user_templates.json", JSON.stringify(window._userTemplates))
    }

    // Mirrors externally-pushed PropertyPanel values onto the floating
    // instance too. Canvas-derived display values (most of PropertyPanel's
    // own computed bindings) don't need this -- both instances already
    // share the same `canvas: activeCanvas` binding and update themselves.
    // This is only for values a signal handler pushes directly onto
    // `propPanel.X`, which the floating instance would otherwise never see.
    function syncFloatingPanel(props) {
        if (!floatingPropPanel) return
        for (const key in props) floatingPropPanel[key] = props[key]
    }

    Settings {
        id: uiSettings
        category: "ui"
        property bool compactDensity: true
        property bool cpkColorsEnabled: true
        property string recentFilesJoined: ""
        property bool fgFullStructure: true
        property string panelState: "docked"
        property real floatingPanelX: -1
        property real floatingPanelY: -1
    }
    Component.onCompleted: {
        Theme.compactDensity = uiSettings.compactDensity
        Theme.cpkColorsEnabled = uiSettings.cpkColorsEnabled
        AppController.fgFullStructure = uiSettings.fgFullStructure
        
        loadUserTemplates()
        
        const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
        fileIO.mkpath(appData + "/recovery")
        // Look for both new .qtsketch format and legacy .ket format
        var recoveredQtsketch = fileIO.listFiles(appData + "/recovery", "autosave_*.qtsketch")
        var recoveredKet = fileIO.listFiles(appData + "/recovery", "autosave_*.ket")
        var recovered = recoveredQtsketch.concat(recoveredKet)
        if (recovered && recovered.length > 0) {
            window._recoveredFiles = recovered
            // Deferred, not opened inline: MessageDialog is a real top-level
            // platform window that centers itself on its parent's geometry. At
            // Component.onCompleted the window has not been shown or sized yet,
            // so it centered on (0,0) and opened at (-208,-137)-(210,110) --
            // title bar and Yes button off-screen, only "No" reachable, while
            // still holding modal keyboard focus (verified via the live window
            // rect). Waiting for the window to actually be sized fixes it.
            recoveryPromptTimer.start()
        }
    }

    Timer {
        id: recoveryPromptTimer
        interval: 250
        repeat: false
        onTriggered: messageDialogsGroup.recoveryDialog.open()
    }

    // Scrolls the workspace so the page's center is centered in the visible
    // viewport. ChemCanvas.fitToMolecule() centers a molecule within the whole
    // A4 page (794x1123 at 100%), but the page is taller than the viewport, so
    // at the default scroll position a perfectly page-centered structure still
    // renders near the bottom edge -- which is what "SMILES import drops the
    // structure off the page" actually was. Deferred by a frame because the
    // content geometry is not final in the same tick a load completes.
    function centerViewOnPage() {
        Qt.callLater(function() {
            const flick = scrollView.contentItem
            if (!flick) return
            // Scroll extent comes from workspaceContent, not flick.contentHeight:
            // ScrollView sizes itself from its content item's implicit geometry and
            // leaves the Flickable's own contentWidth/contentHeight at -1, which
            // silently clamped every target scroll position to 0.
            const maxY = Math.max(0, workspaceContent.height - flick.height)
            const maxX = Math.max(0, workspaceContent.width  - flick.width)
            const targetY = docRect.y + (docRect.height * window.zoomLevel) / 2 - flick.height / 2
            const targetX = docRect.x + (docRect.width  * window.zoomLevel) / 2 - flick.width  / 2
            flick.contentY = Math.max(0, Math.min(targetY, maxY))
            flick.contentX = Math.max(0, Math.min(targetX, maxX))
        })
    }

    Connections {
        target: AppController
        function onFgFullStructureChanged() { uiSettings.fgFullStructure = AppController.fgFullStructure }
    }

    visible: true
    width: 1400
    height: 900
    title: "sketch" + (activeCanvas && activeCanvas.isDirty ? " *" : "")
    color: Theme.background

    property bool _forceQuit: false

    onClosing: function(close_event) {
        if (window._forceQuit) return
        // Check ALL tabs for unsaved changes, not just the active one
        var anyDirty = false
        for (var i = 0; i < canvasRepeater.count; i++) {
            var item = canvasRepeater.itemAt(i)
            if (item && item.isDirty) { anyDirty = true; break }
        }
        if (anyDirty) {
            close_event.accepted = false
            messageDialogsGroup.unsavedChangesDialog.open()
        } else {
            // Clean exit
            const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
            // Clean up both new and old recovery formats
            var recoveredQtsketch = fileIO.listFiles(appData + "/recovery", "autosave_*.qtsketch")
            var recoveredKet = fileIO.listFiles(appData + "/recovery", "autosave_*.ket")
            var allRecovered = recoveredQtsketch.concat(recoveredKet)
            for (var j = 0; j < allRecovered.length; j++) {
                fileIO.remove(appData + "/recovery/" + allRecovered[j])
            }
        }
    }

    // Change cursor globally when processing
    ProgressBar {
        id: processingBar
        anchors.top: parent.top
        anchors.left: parent.left
        anchors.right: parent.right
        height: 2
        indeterminate: true
        visible: window.isProcessing
        z: 9999
    }

    ColumnLayout {
        anchors.fill: parent
        spacing: 0

        // Indigo Backend
        IndigoService {
            id: indigoSvc
        }

        // Imago OCR Backend — chemical-structure image recognition for Insert Image
        ImagoService {
            id: imagoSvc
            onImageRecognized: function(molfile, warningsCount, error) {
                window.isProcessing = false
                if (!error && molfile && warningsCount <= fileDialogsGroup.imageFileDialog.maxAcceptableWarnings && activeSketch) {
                    activeSketch.sendCommand("insertRecognizedStructure", [molfile, fileDialogsGroup.imageFileDialog.chemX, fileDialogsGroup.imageFileDialog.chemY])
                    imageConfirmBanner.showBanner(warningsCount)
                } else {
                    // Recognition failed or low-confidence — fall back to a plain image embed.
                    const dataUri = fileIO.readImageAsDataUri(fileDialogsGroup.imageFileDialog.pendingFileUrl)
                    if (dataUri === "") {
                        messageDialogsGroup.workerErrorDialog.errorText = "Failed to load image. Ensure it is a supported format (png, jpg, gif, bmp) and under 5 MB."
                        messageDialogsGroup.workerErrorDialog.severe = false
                        messageDialogsGroup.workerErrorDialog.open()
                    } else if (activeSketch) {
                        activeSketch.addImage(dataUri, fileDialogsGroup.imageFileDialog.chemX, fileDialogsGroup.imageFileDialog.chemY, 1.5, 1.5)
                    }
                }
            }
        }

        // Global Shortcuts
        Shortcut { sequence: "S"; onActivated: if (activeCanvas) activeCanvas.currentTool = "SELECT" }
        Shortcut { sequence: "E"; onActivated: if (activeCanvas) activeCanvas.currentTool = "ERASE" }
        Shortcut { sequence: "B"; onActivated: if (activeCanvas) activeCanvas.currentTool = "BOND_1" }
        Shortcut { sequence: "R"; onActivated: if (activeCanvas) activeCanvas.currentTool = "TEMPLATE_BENZENE" }
        Shortcut { sequence: "Ctrl+R"; onActivated: rgroupPanel.open() }
        Shortcut { sequence: "Ctrl+K"; onActivated: commandPalette.open() }

        CommandPalette {
            id: commandPalette
            window: window
        }

        Connections {
            target: indigoSvc
            function onLayoutFinished(newMol, error) {
                window.isProcessing = false
                // IndigoService::layout preserves the structure's original centroid,
                // so no forced re-centering here -- don't disturb where the user had it.
                if (newMol && activeCanvas) {
                    activeCanvas.loadMolfile(newMol)
                } else if (error) {
                    messageDialogsGroup.workerErrorDialog.errorText = "Layout: " + error
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onCommonScaffoldFinished(result, error) {
                window.isProcessing = false
                if (result && activeCanvas) {
                    activeCanvas.loadMolfile(result)
                } else {
                    messageDialogsGroup.workerErrorDialog.errorText = "Find Common Scaffold: " + (error || "unknown error")
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onRgroupDecompositionFinished(result, error) {
                window.isProcessing = false
                if (result && activeCanvas) {
                    activeCanvas.loadMolfile(result)
                } else {
                    messageDialogsGroup.workerErrorDialog.errorText = "Decompose to R-Groups: " + (error || "unknown error")
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onRgroupPerMoleculeDecompositionFinished(resultsJson, error) {
                window.isProcessing = false
                if (resultsJson) {
                    window._rgroupPerMoleculeResults = JSON.parse(resultsJson)
                    taskDialogsGroup.rgroupPerMoleculeResultsDialog.open()
                } else {
                    messageDialogsGroup.workerErrorDialog.errorText = "Decompose to R-Groups (per-molecule): " + (error || "unknown error")
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onSimilarityRankFinished(resultJson, error) {
                window.isProcessing = false
                if (resultJson) {
                    try {
                        const arr = JSON.parse(resultJson)
                        arr.sort((a, b) => b.score - a.score)
                        const labels = window._pendingBatchLabels || []
                        const lines = arr.map((r, i) => {
                            const name = labels[r.index] || ("Record " + (r.index + 1))
                            return (i + 1) + ". " + name + " — " + (r.score * 100).toFixed(1) + "%"
                        })
                        messageDialogsGroup.similarityRankResultDialog.text = lines.join("\n")
                        messageDialogsGroup.similarityRankResultDialog.open()
                    } catch (e) {
                        messageDialogsGroup.workerErrorDialog.errorText = "Rank by Similarity: failed to parse results."
                        messageDialogsGroup.workerErrorDialog.severe = false
                        messageDialogsGroup.workerErrorDialog.open()
                    }
                } else {
                    messageDialogsGroup.workerErrorDialog.errorText = "Rank by Similarity: " + (error || "unknown error")
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onBatchAlignFinished(resultJson, error) {
                window.isProcessing = false
                if (resultJson && activeSketch) {
                    activeSketch.sendCommand("realignSdfBatch", [resultJson])
                } else {
                    messageDialogsGroup.workerErrorDialog.errorText = "Align Batch to Common Scaffold: " + (error || "unknown error")
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onRdfBatchParsed(recordsJson, error) {
                if (recordsJson && activeSketch) {
                    activeSketch.sendCommand("deserializeRdfBatch", [recordsJson])
                } else {
                    window.isProcessing = false
                    messageDialogsGroup.workerErrorDialog.errorText = "Open RDF: " + (error || "unknown error")
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onIndigoBatchParsed(recordsJson, error) {
                if (recordsJson && activeSketch) {
                    activeSketch.sendCommand("deserializeIndigoBatch", [recordsJson])
                } else {
                    window.isProcessing = false
                    messageDialogsGroup.workerErrorDialog.errorText = "Open batch file: " + (error || "unknown error")
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onClean2dFinished(newMol) {
                window.isProcessing = false
                if (newMol && activeCanvas) activeCanvas.loadMolfile(newMol)
            }
            function onAromatizeFinished(newMol) {
                window.isProcessing = false
                if (newMol && activeCanvas) activeCanvas.loadMolfile(newMol)
            }
            function onDearomatizeFinished(newMol) {
                window.isProcessing = false
                if (newMol && activeCanvas) activeCanvas.loadMolfile(newMol)
            }
            function onUnfoldHydrogensFinished(newMol) {
                window.isProcessing = false
                if (newMol && activeCanvas) activeCanvas.loadMolfile(newMol)
            }
            function onFoldHydrogensFinished(newMol) {
                window.isProcessing = false
                if (newMol && activeCanvas) activeCanvas.loadMolfile(newMol)
            }
            function onNormalizeFinished(newMol) {
                window.isProcessing = false
                if (newMol && activeCanvas) activeCanvas.loadMolfile(newMol)
            }
            function onStandardizeFinished(newMol) {
                window.isProcessing = false
                if (newMol && activeCanvas) activeCanvas.loadMolfile(newMol)
            }
            function onSmilesFinished(smiles) {
                if (smiles) indigoSvc.copyToClipboard(smiles)
            }
            function onCanonicalSmilesFinished(smiles) {
                if (smiles) indigoSvc.copyToClipboard(smiles)
            }
            function onInchiFinished(inchi) {
                if (inchi) indigoSvc.copyToClipboard(inchi)
            }
            function onInchiKeyFinished(inchiKey) {
                if (inchiKey) indigoSvc.copyToClipboard(inchiKey)
            }
            function onHashFinished(hash) {
                if (hash) indigoSvc.copyToClipboard(hash)
            }
            function onSimilarityFinished(result) {
                if (result) {
                    messageDialogsGroup.similarityResultDialog.text = "Tanimoto similarity: " + result
                    messageDialogsGroup.similarityResultDialog.open()
                }
            }
            function onMassCompositionFinished(result) {
                if (result) indigoSvc.copyToClipboard(result)
            }
            function onPkaValuesFinished(result) {
                if (result) indigoSvc.copyToClipboard(result)
            }
            function onRenderFinished(success, error) {
                if (!success) {
                    const isPdf = window.pendingRenderUrl.toString().toLowerCase().endsWith(".pdf")
                    messageDialogsGroup.workerErrorDialog.errorText = (isPdf ? "PDF export failed: " : "SVG export failed: ") + error
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onReactionMappingFinished(result, error) {
                window.isProcessing = false
                if (result && activeCanvas) {
                    activeCanvas.loadMolfile(result)
                } else {
                    messageDialogsGroup.workerErrorDialog.errorText = "Atom Mapping: " + (error || "unknown error")
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onIonizeFinished(result, error) {
                window.isProcessing = false
                if (result && activeCanvas) {
                    activeCanvas.loadMolfile(result)
                } else {
                    messageDialogsGroup.workerErrorDialog.errorText = "Ionize at pH: " + (error || "unknown error")
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            function onPropertiesReady(mw, mono, mf, atoms, bonds, tpsa, logp, hba, hbd, rotBonds, molarRefractivity, pka, heavyAtoms, isChiral, mostAbundantMass, fragmentCount, ringCount) {
                propPanel.molMW       = mw
                propPanel.molMono     = mono
                propPanel.molFormula  = mf
                propPanel.molAtoms    = atoms
                propPanel.molBonds    = bonds
                propPanel.molTPSA     = tpsa
                propPanel.molLogP     = logp
                propPanel.molHBA      = hba
                propPanel.molHBD      = hbd
                propPanel.molRotBonds = rotBonds
                propPanel.molMolarRefractivity = molarRefractivity
                propPanel.molPka      = pka
                propPanel.molHeavyAtoms = heavyAtoms
                propPanel.molIsChiral   = isChiral
                propPanel.molMostAbundantMass = mostAbundantMass
                propPanel.molFragmentCount = fragmentCount
                propPanel.molRingCount = ringCount
                syncFloatingPanel({
                    molMW: mw, molMono: mono, molFormula: mf, molAtoms: atoms,
                    molBonds: bonds, molTPSA: tpsa, molLogP: logp, molHBA: hba,
                    molHBD: hbd, molRotBonds: rotBonds, molMolarRefractivity: molarRefractivity,
                    molPka: pka, molHeavyAtoms: heavyAtoms, molIsChiral: isChiral,
                    molMostAbundantMass: mostAbundantMass, molFragmentCount: fragmentCount,
                    molRingCount: ringCount
                })
            }
            function onStereoDescriptorsReady(jsonMap) {
                if (activeSketch) activeSketch.setStereoDescriptors(jsonMap)
            }
            function onCheckIssuesReady(structuredJson) {
                if (activeSketch) activeSketch.setCheckIssues(structuredJson)
            }
            function onSubstructureSearchFinished(resultJson) {
                const parsed = JSON.parse(resultJson)
                if (parsed.error) {
                    smartsSearchPopup.statusIsError = true
                    smartsSearchPopup.statusText = parsed.error
                    return
                }
                smartsSearchPopup.statusIsError = false
                smartsSearchPopup.statusText = parsed.matchCount + " match" + (parsed.matchCount === 1 ? "" : "es") + (parsed.truncated ? " (showing first 200)" : "") + " found"
                if (activeSketch) activeSketch.sendCommand("selectSubstructureMatches", [JSON.stringify(parsed)])
            }
            function onCheckFinished(report) {
                // Inline check warnings (badges, red bonds) are driven by checkIssuesReady,
                // which fires silently on every calc_props cycle. The modal dialog is
                // explicit-only: opened solely when the user clicks the Validate button.
                if (window.explicitCheckPending) {
                    window.explicitCheckPending = false
                    taskDialogsGroup.checkResultDialog.reportText = report
                    taskDialogsGroup.checkResultDialog.open()
                }
            }
            function onIupacNameReady(name, error) {
                propPanel.iupacLoading = false
                propPanel.iupacName = name
                propPanel.iupacError = error
                syncFloatingPanel({ iupacLoading: false, iupacName: name, iupacError: error })
            }

            // Biopolymer load
            function onBiopolymerLoaded(molfile) {
                window.isProcessing = false
                if (activeCanvas) activeCanvas.loadMolfile(molfile)
                biopolymerDialog.onLoadSuccess()
            }
            function onBiopolymerLoadError(error) {
                window.isProcessing = false
                biopolymerDialog.showLoadError(error)
            }

            // Biopolymer export results → shown in dialog textarea
            function onBioSequenceReady(text) { biopolymerDialog.showExportResult("Sequence", text) }
            function onBioFastaReady(text)    { biopolymerDialog.showExportResult("FASTA",    text) }
            function onBioHelmReady(text)     { biopolymerDialog.showExportResult("HELM",     text) }
            function onBioIdtReady(text)      { biopolymerDialog.showExportResult("IDT",      text) }
            function onBioAxoLabsReady(text)  { biopolymerDialog.showExportResult("AxoLabs",  text) }
        }

        Timer {
            id: propUpdateTimer
            interval: 600
            repeat: false
            onTriggered: if (activeSketch) activeSketch.requestSerialize("calc_props")
        }

        Timer {
            id: autosaveTimer
            interval: 60000
            repeat: true
            running: true
            onTriggered: {
                for (let i = 0; i < DocumentManager.docIds.length; i++) {
                    let docId = DocumentManager.docIds[i]
                    if (window.dirtyDocs[docId]) {
                        const sketch = DocumentManager.documentFor(docId)
                        if (sketch) {
                            // Use new project serialization format instead of KET
                            sketch.requestProjectSerialization("autosave")
                        }
                    }
                }
            }
        }

        Connections {
            target: window.activeSketch
            readonly property var structureResponseRoutes: ({
                save: (data) => fileIO.write(window.pendingSaveUrl, data),
                autosave: (data) => {
                    const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
                    fileIO.write(appData + "/recovery/autosave_" + window.activeSketch.modelData + ".qtsketch", data)
                },
                render_svg: (data) => indigoSvc.renderToFile(data, window.pendingRenderUrl, "svg"),
                render_pdf: (data) => indigoSvc.renderToFile(data, window.pendingRenderUrl, "pdf"),
                render_png: (data) => indigoSvc.renderToFile(data, window.pendingRenderUrl, "png"),
                render_grid: (data) => indigoSvc.renderReactionGridToFile(data, window.pendingRenderUrl, "pdf"),
                automap: (data) => indigoSvc.autoMapReaction(data),
                ionize: (data) => indigoSvc.ionizeAtPh(data, window.pendingIonizePh),
                clear_mapping: (data) => indigoSvc.clearReactionMapping(data),
                correct_reacting_centers: (data) => indigoSvc.correctReactingCenters(data),
                smiles: (data) => indigoSvc.smiles(data),
                canonical_smiles: (data) => indigoSvc.canonicalSmiles(data),
                inchi: (data) => indigoSvc.inchi(data),
                inchikey: (data) => indigoSvc.inchiKey(data),
                hash: (data) => indigoSvc.hash(data),
                similarity: (data) => indigoSvc.similarity(data, window.pendingSimilarityRef),
                mass_composition: (data) => indigoSvc.massComposition(data),
                pka_values: (data) => indigoSvc.pkaValues(data),
                layout: (data) => { window.isProcessing = true; indigoSvc.layout(data); },
                clean2d: (data) => { window.isProcessing = true; indigoSvc.clean2d(data); },
                aromatize: (data) => { window.isProcessing = true; indigoSvc.aromatize(data); },
                dearomatize: (data) => { window.isProcessing = true; indigoSvc.dearomatize(data); },
                unfoldH: (data) => { window.isProcessing = true; indigoSvc.unfoldHydrogens(data); },
                foldH: (data) => { window.isProcessing = true; indigoSvc.foldHydrogens(data); },
                normalize: (data) => { window.isProcessing = true; indigoSvc.normalize(data); },
                standardize: (data) => { window.isProcessing = true; indigoSvc.standardize(data); },
                check: (data) => indigoSvc.checkStructure(data),
                iupac_name: (data) => indigoSvc.generateIupacName(data),
                calc_props: (data) => {
                    indigoSvc.calcProperties(data);
                    indigoSvc.calcStereoDescriptors(data);
                    indigoSvc.checkStructure(data);
                    if (activeSketch) {
                        activeSketch.sendCommand("getMoleculeName", []);
                        activeSketch.sendCommand("getSdfProps", []);
                    }
                },
                bio_seq: (data) => indigoSvc.exportBioSequence(data),
                bio_fasta: (data) => indigoSvc.exportBioFasta(data),
                bio_helm: (data) => indigoSvc.exportBioHelm(data),
                bio_idt: (data) => indigoSvc.exportBioIdt(data),
                bio_axolabs: (data) => indigoSvc.exportBioAxoLabs(data),
                atom_props: (data) => {
                    const props = JSON.parse(data);
                    if (props && props.id !== undefined)
                        atomPropsDialog.openForAtom(props.id, props);
                },
                clipboard_ket: (data) => {
                    if (data && data.length > 0 && activeSketch)
                        activeSketch.setOsClipboardText(data);
                },
                sdf_batch_list: (data) => {
                    window.isProcessing = false;
                    const parsed = JSON.parse(data);
                    if (parsed.count <= 1) {
                        if (activeSketch) activeSketch.sendCommand("loadSdfBatchRecord", [0]);
                    } else {
                        sdfRecordPicker.recordCount = parsed.count;
                        sdfRecordPicker.model = parsed.records;
                        sdfRecordPicker.open();
                    }
                },
                sdf_batch_realigned: (data) => {
                    const parsed = JSON.parse(data);
                    sdfRecordPicker.recordCount = parsed.count;
                    sdfRecordPicker.model = parsed.records;
                },
                sdf_batch_load: (data) => {
                    window.isProcessing = false;
                    if (activeCanvas) {
                        activeCanvas._needsCentering = true;
                        activeCanvas.refresh();
                        window.centerViewOnPage();
                        setDocFile(DocumentManager.activeDocId, pendingSdfBatchUrl);
                    }
                },
                sdf_batch_molfiles: (data) => {
                    const parsed = JSON.parse(data);
                    if (parsed.molfiles && parsed.molfiles.length >= 2) {
                        if (window._pendingBatchAction === "export_grid") {
                            window._pendingBatchGridMolfiles = parsed.molfiles;
                            fileDialogsGroup.batchGridSaveDialog.open();
                        } else if (window._pendingBatchAction === "export_file") {
                            window._pendingBatchGridMolfiles = parsed.molfiles;
                            fileDialogsGroup.batchFileSaveDialog.open();
                        } else if (window._pendingBatchAction === "align") {
                            indigoSvc.alignBatchToScaffold(parsed.molfiles);
                        } else if (window._pendingBatchAction === "decompose") {
                            indigoSvc.decomposeToRGroups(parsed.molfiles);
                        } else if (window._pendingBatchAction === "decompose_permolecule") {
                            indigoSvc.decomposeToRGroupsPerMolecule(parsed.molfiles, parsed.labels);
                        } else if (window._pendingBatchAction === "similarity") {
                            window._pendingBatchLabels = parsed.labels || [];
                            if (!window.pendingSimilarityRefMolfile) {
                                window.isProcessing = false;
                                messageDialogsGroup.workerErrorDialog.errorText = "Rank by Similarity: no active structure to compare against.";
                                messageDialogsGroup.workerErrorDialog.severe = false;
                                messageDialogsGroup.workerErrorDialog.open();
                            } else {
                                indigoSvc.rankBySimilarity(window.pendingSimilarityRefMolfile, parsed.molfiles);
                            }
                        } else {
                            indigoSvc.findCommonScaffold(parsed.molfiles);
                        }
                    } else {
                        window.isProcessing = false;
                        messageDialogsGroup.workerErrorDialog.errorText = "Batch structure analysis: not enough valid structures in this batch.";
                        messageDialogsGroup.workerErrorDialog.severe = false;
                        messageDialogsGroup.workerErrorDialog.open();
                    }
                },
                smarts_search: (data) => {
                    indigoSvc.substructureSearch(data, window.pendingSmartsQuery);
                },
                mol_name: (data) => {
                    propPanel.molName = data;
                    syncFloatingPanel({ molName: data })
                },
                sdf_props: (data) => {
                    propPanel.sdfProps = JSON.parse(data);
                    syncFloatingPanel({ sdfProps: JSON.parse(data) })
                },
                clipboard_preview: (data) => {
                    if (activeCanvas) activeCanvas.clipboardPreview = data === "null" ? null : JSON.parse(data);
                },
                biopolymer_seq_view: (data) => {
                    biopolymerSeqView.openSnapshot(data);
                }
            })

            function onStructureReady(reqId, data) {
                const handler = structureResponseRoutes[reqId];
                if (handler) handler(data);
            }
            function onStateUpdated(state, selection, dirty, undoState, redoState, result) {
                if (result !== "stereoUpdated" && result !== "checkUpdated" && result !== "selectionOnly") {
                    propUpdateTimer.restart()
                    // A generated name is only valid for the exact structure it was
                    // generated from; any real structural edit (or undo/redo) invalidates it.
                    // Pure selection changes (click/rect/lasso select, select-all, etc.) are
                    // tagged "selectionOnly" by V8Process and don't touch the structure, so
                    // they must not wipe a name the user just generated.
                    propPanel.iupacName = ""
                    propPanel.iupacError = ""
                    syncFloatingPanel({ iupacName: "", iupacError: "" })
                }
            }
            function onErrorOccurred(error) {
                window.isProcessing = false
                propPanel.iupacLoading = false
                syncFloatingPanel({ iupacLoading: false })
                messageDialogsGroup.workerErrorDialog.errorText = error
                // Explicit even though severe defaults to true -- guards against
                // a previous recoverable-error call site's severe = false
                // lingering on this shared dialog instance.
                messageDialogsGroup.workerErrorDialog.severe = true
                messageDialogsGroup.workerErrorDialog.open()
            }
        }

        AppMenus {
            win: window
            Layout.fillWidth: true
        }

        // Top Toolbar
        MainToolbar {
            id: mainToolbar
            win: window
        }

        // Horizontal separator
        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 1
            color: Theme.outline
        }

        // Document tabs
        RowLayout {
            Layout.fillWidth: true
            spacing: 4

            Item {
                id: docTabBarWrapper
                Layout.fillWidth: true
                implicitHeight: docTabBar.implicitHeight

                TabBar {
                    id: docTabBar
                    anchors.fill: parent
                    currentIndex: DocumentManager.docIds.indexOf(DocumentManager.activeDocId)
                    onCurrentIndexChanged: {
                        var ids = DocumentManager.docIds
                        if (currentIndex >= 0 && currentIndex < ids.length) DocumentManager.activeDocId = ids[currentIndex]
                    }

                    Repeater {
                        model: DocumentManager.docIds
                        TabButton {
                            id: tabBtn
                            indicator: Item {}
                            required property int modelData
                            // titleRev forces re-evaluation after open / save-as / rename
                            text: (window.titleRev, (window.dirtyDocs[modelData] ? "● " : "") + window.titleFor(modelData))

                            // Content-sized and left-aligned, the Fluent tab behaviour.
                            // TabBar stretches its buttons to fill the bar by default,
                            // so a lone document tab spanned the whole window with its
                            // label centred in the middle of empty space. The +32 is
                            // room for the close affordance overlaid on the right.
                            width: Math.min(240, Math.max(120, implicitWidth + 32))

                            onDoubleClicked: {
                                taskDialogsGroup.renameDialog.docId = modelData
                                taskDialogsGroup.renameDialog.open()
                            }

                            Text {
                                // Segoe Fluent Icons Dismiss (U+E8BB) -- same glyph the
                                // window's own Close caption button uses.
                                text: String.fromCharCode(0xE8BB)
                                anchors.right: parent ? parent.right : undefined
                                anchors.verticalCenter: parent ? parent.verticalCenter : undefined
                                anchors.rightMargin: 8
                                font { family: Theme.fontIcons; pixelSize: 10 }
                                color: tabMouseArea.containsMouse ? Theme.error : Theme.textSecondary
                                MouseArea {
                                    id: tabMouseArea
                                    anchors.fill: parent
                                    hoverEnabled: true
                                    onClicked: {
                                        if (DocumentManager.docIds.length <= 1) return
                                        if (window.dirtyDocs[tabBtn.modelData]) {
                                            // Switch to the tab being closed first: fileDialogsGroup.saveDialog/the
                                            // structureReady Connections below only ever operate on
                                            // activeSketch/activeCanvas, so saving a *different*,
                                            // still-background tab here would silently save the
                                            // wrong document's content instead.
                                            DocumentManager.activeDocId = tabBtn.modelData
                                            window._pendingCloseDocId = tabBtn.modelData
                                            messageDialogsGroup.unsavedChangesDialog.open()
                                        } else {
                                            DocumentManager.closeDocument(tabBtn.modelData)
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                Rectangle {
                    id: tabIndicator
                    height: 2
                    color: Theme.accent
                    y: docTabBar.height - height
                    x: docTabBar.currentItem ? docTabBar.currentItem.x : 0
                    width: docTabBar.currentItem ? docTabBar.currentItem.width : 0
                    Behavior on x { NumberAnimation { duration: Theme.durationFast; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }
                    Behavior on width { NumberAnimation { duration: Theme.durationFast; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }
                }
            }

            Button {
                text: "+"
                onClicked: DocumentManager.addDocument()
            }
        }

        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 0

            ToolPanel {
                Layout.preferredWidth: Theme.toolPanelWidth
                Layout.fillHeight: true
                sketch: window.activeSketch
                win: window
                currentTool: activeCanvas ? activeCanvas.currentTool : "SELECT"
                onToolSelected: (toolId) => {
                    if (!activeCanvas) return
                    if (toolId === "ATOM_ANY") {
                        periodicTablePopup.open()
                    } else {
                        if (activeCanvas.currentTool === toolId) {
                            activeCanvas.currentTool = "SELECT"
                        } else {
                            activeCanvas.currentTool = toolId
                        }
                    }
                }
            }

            // Thin separator
            Rectangle {
                Layout.preferredWidth: 1
                Layout.fillHeight: true
                color: Theme.outline
            }

            ColumnLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                spacing: 0

                // Top Ruler
                Canvas {
                    id: topRuler
                    antialiasing: true
                    Layout.fillWidth: true
                    Layout.preferredHeight: 24
                    Layout.leftMargin: 24 // Align with the left ruler width

                    // Ruler units label, drawn in the empty corner this leftMargin
                    // leaves above leftRuler -- topRuler doesn't clip its children,
                    // so a negative x here still renders inside that gap.
                    Text {
                        text: "cm"
                        x: -20
                        y: 6
                        color: Theme.rulerColor
                        font { pixelSize: Theme.fontSizeCaption; family: Theme.fontMono }
                    }

                    property real scrollX: scrollView.contentItem.contentX
                    property real zoom: window.zoomLevel
                    property real docX: docRect.x
                    property var margins: window.pageMargins
                    property string pageSize: mainToolbar.pageSizeCombo.currentText
                    
                    readonly property real cmPixels: window.pixelsPerCm * zoom
                    readonly property real startX: docX - scrollX
                    
                    onScrollXChanged: requestPaint()
                    onZoomChanged: requestPaint()
                    onDocXChanged: requestPaint()
                    onWidthChanged: requestPaint()
                    onMarginsChanged: requestPaint()
                    onPageSizeChanged: requestPaint()
                    
                    onPaint: {
                        const ctx = getContext("2d");
                        ctx.clearRect(0, 0, width, height);
                        ctx.fillStyle = Theme.surface.toString();
                        ctx.fillRect(0, 0, width, height);
                        ctx.strokeStyle = Theme.rulerColor;
                        ctx.fillStyle = Theme.rulerColor;
                        ctx.font = Theme.fontSizeCaption + "px " + Theme.fontMonoCss;
                        ctx.beginPath();
                        
                        for (let i = -50; i < 150; i += 0.1) {
                            const x = startX + (i * cmPixels);
                            if (x < 0) continue;
                            if (x > width) break;
                            
                            if (Math.abs(i - Math.round(i)) < 0.01) {
                                ctx.moveTo(x, 0); ctx.lineTo(x, height);
                                if (i !== 0) ctx.fillText(Math.round(i).toString(), x + 2, 10);
                            } else if (Math.abs(i - Math.round(i*2)/2) < 0.01) {
                                ctx.moveTo(x, height/2); ctx.lineTo(x, height);
                            } else {
                                ctx.moveTo(x, height - 5); ctx.lineTo(x, height);
                            }
                        }
                        ctx.stroke();
                        
                        // Shade left and right margins
                        ctx.save();
                        ctx.globalAlpha = 0.15;
                        ctx.fillStyle = Theme.rulerColor;
                        
                        const pWidthMm = window.pageSizeMm(mainToolbar.pageSizeCombo.currentText).w;
                        const marginLeftPx = startX + (window.pageMargins.left * cmPixels);
                        const marginRightPx = startX + ((pWidthMm/10 - window.pageMargins.right) * cmPixels);
                        const pageWidthPx = startX + ((pWidthMm/10) * cmPixels);
                        
                        const leftRectW = Math.max(0, Math.min(marginLeftPx, width));
                        ctx.fillRect(0, 0, leftRectW, height);
                        
                        const rightRectX = Math.max(0, Math.min(marginRightPx, width));
                        const rightRectW = Math.max(0, Math.min(pageWidthPx, width) - rightRectX);
                        ctx.fillRect(rightRectX, 0, rightRectW, height);
                        
                        ctx.restore();
                        
                        ctx.beginPath();
                        ctx.moveTo(0, height);
                        ctx.lineTo(width, height);
                        ctx.stroke();
                    }

                    MouseArea {
                        id: leftMarginDrag
                        width: 6
                        height: parent.height
                        x: parent.startX + (window.pageMargins.left * parent.cmPixels) - 3
                        y: 0
                        cursorShape: Qt.SizeHorCursor
                        
                        property real dragStartX: 0
                        property real marginStartX: 0
                        
                        onPressed: (mouse) => {
                            let parentPos = mapToItem(parent, mouse.x, mouse.y)
                            dragStartX = parentPos.x
                            marginStartX = window.pageMargins.left
                        }
                        
                        onPositionChanged: (mouse) => {
                            if (pressed) {
                                let parentPos = mapToItem(parent, mouse.x, mouse.y)
                                let deltaX = parentPos.x - dragStartX
                                let deltaCm = deltaX / parent.cmPixels
                                let newMargin = marginStartX + deltaCm
                                
                                let pWidthCm = window.pageSizeMm(mainToolbar.pageSizeCombo.currentText).w / 10
                                let maxMargin = pWidthCm - window.pageMargins.right - 1.0
                                newMargin = Math.max(0.0, Math.min(newMargin, maxMargin))
                                
                                let margins = {
                                    left: newMargin,
                                    right: window.pageMargins.right,
                                    top: window.pageMargins.top,
                                    bottom: window.pageMargins.bottom
                                }
                                window.pageMargins = margins
                            }
                        }
                    }

                    MouseArea {
                        id: rightMarginDrag
                        width: 6
                        height: parent.height
                        x: {
                            let pWidthCm = window.pageSizeMm(mainToolbar.pageSizeCombo.currentText).w / 10
                            return parent.startX + ((pWidthCm - window.pageMargins.right) * parent.cmPixels) - 3
                        }
                        y: 0
                        cursorShape: Qt.SizeHorCursor
                        
                        property real dragStartX: 0
                        property real marginStartX: 0
                        
                        onPressed: (mouse) => {
                            let parentPos = mapToItem(parent, mouse.x, mouse.y)
                            dragStartX = parentPos.x
                            marginStartX = window.pageMargins.right
                        }
                        
                        onPositionChanged: (mouse) => {
                            if (pressed) {
                                let parentPos = mapToItem(parent, mouse.x, mouse.y)
                                let deltaX = parentPos.x - dragStartX
                                let deltaCm = deltaX / parent.cmPixels
                                let newMargin = marginStartX - deltaCm
                                
                                let pWidthCm = window.pageSizeMm(mainToolbar.pageSizeCombo.currentText).w / 10
                                let maxMargin = pWidthCm - window.pageMargins.left - 1.0
                                newMargin = Math.max(0.0, Math.min(newMargin, maxMargin))
                                
                                let margins = {
                                    left: window.pageMargins.left,
                                    right: newMargin,
                                    top: window.pageMargins.top,
                                    bottom: window.pageMargins.bottom
                                }
                                window.pageMargins = margins
                            }
                        }
                    }
                }

                RowLayout {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    spacing: 0

                    // Left Ruler
                    Canvas {
                        id: leftRuler
                        antialiasing: true
                        Layout.preferredWidth: 24
                        Layout.fillHeight: true
                        
                        property real scrollY: scrollView.contentItem.contentY
                        property real zoom: window.zoomLevel
                        property real docY: docRect.y
                        property var margins: window.pageMargins
                        property string pageSize: mainToolbar.pageSizeCombo.currentText
                        
                        readonly property real cmPixels: window.pixelsPerCm * zoom
                        readonly property real startY: docY - scrollY
                        
                        onScrollYChanged: requestPaint()
                        onZoomChanged: requestPaint()
                        onDocYChanged: requestPaint()
                        onHeightChanged: requestPaint()
                        onMarginsChanged: requestPaint()
                        onPageSizeChanged: requestPaint()
                        
                        onPaint: {
                            const ctx = getContext("2d");
                            ctx.clearRect(0, 0, width, height);
                            ctx.fillStyle = Theme.surface.toString();
                            ctx.fillRect(0, 0, width, height);
                            ctx.strokeStyle = Theme.rulerColor;
                            ctx.fillStyle = Theme.rulerColor;
                            ctx.font = Theme.fontSizeCaption + "px " + Theme.fontMonoCss;
                            ctx.beginPath();
                            
                            for (let i = -50; i < 150; i += 0.1) {
                                const y = startY + (i * cmPixels);
                                if (y < 0) continue;
                                if (y > height) break;
                                
                                if (Math.abs(i - Math.round(i)) < 0.01) {
                                    ctx.moveTo(0, y); ctx.lineTo(width, y);
                                    if (i !== 0) ctx.fillText(Math.round(i).toString(), 2, y + 10);
                                } else if (Math.abs(i - Math.round(i*2)/2) < 0.01) {
                                    ctx.moveTo(width/2, y); ctx.lineTo(width, y);
                                } else {
                                    ctx.moveTo(width - 5, y); ctx.lineTo(width, y);
                                }
                            }
                            ctx.stroke();
                            
                            // Shade top and bottom margins
                            ctx.save();
                            ctx.globalAlpha = 0.15;
                            ctx.fillStyle = Theme.rulerColor;
                            
                            const pHeightMm = window.pageSizeMm(mainToolbar.pageSizeCombo.currentText).h;
                            const marginTopPx = startY + (window.pageMargins.top * cmPixels);
                            const marginBottomPx = startY + ((pHeightMm/10 - window.pageMargins.bottom) * cmPixels);
                            const pageHeightPx = startY + ((pHeightMm/10) * cmPixels);
                            
                            const topRectH = Math.max(0, Math.min(marginTopPx, height));
                            ctx.fillRect(0, 0, width, topRectH);
                            
                            const bottomRectY = Math.max(0, Math.min(marginBottomPx, height));
                            const bottomRectH = Math.max(0, Math.min(pageHeightPx, height) - bottomRectY);
                            ctx.fillRect(0, bottomRectY, width, bottomRectH);
                            
                            ctx.restore();
                            
                            ctx.beginPath();
                            ctx.moveTo(width, 0);
                            ctx.lineTo(width, height);
                            ctx.stroke();
                        }

                        MouseArea {
                            id: topMarginDrag
                            width: parent.width
                            height: 6
                            x: 0
                            y: parent.startY + (window.pageMargins.top * parent.cmPixels) - 3
                            cursorShape: Qt.SizeVerCursor
                            
                            property real dragStartY: 0
                            property real marginStartY: 0
                            
                            onPressed: (mouse) => {
                                let parentPos = mapToItem(parent, mouse.x, mouse.y)
                                dragStartY = parentPos.y
                                marginStartY = window.pageMargins.top
                            }
                            
                            onPositionChanged: (mouse) => {
                                if (pressed) {
                                    let parentPos = mapToItem(parent, mouse.x, mouse.y)
                                    let deltaY = parentPos.y - dragStartY
                                    let deltaCm = deltaY / parent.cmPixels
                                    let newMargin = marginStartY + deltaCm
                                    
                                    let pHeightCm = window.pageSizeMm(mainToolbar.pageSizeCombo.currentText).h / 10
                                    let maxMargin = pHeightCm - window.pageMargins.bottom - 1.0
                                    newMargin = Math.max(0.0, Math.min(newMargin, maxMargin))
                                    
                                    let margins = {
                                        left: window.pageMargins.left,
                                        right: window.pageMargins.right,
                                        top: newMargin,
                                        bottom: window.pageMargins.bottom
                                    }
                                    window.pageMargins = margins
                                }
                            }
                        }

                        MouseArea {
                            id: bottomMarginDrag
                            width: parent.width
                            height: 6
                            x: 0
                            y: {
                                const pHeightCm = window.pageSizeMm(mainToolbar.pageSizeCombo.currentText).h / 10
                                return parent.startY + ((pHeightCm - window.pageMargins.bottom) * parent.cmPixels) - 3
                            }
                            cursorShape: Qt.SizeVerCursor
                            
                            property real dragStartY: 0
                            property real marginStartY: 0
                            
                            onPressed: (mouse) => {
                                let parentPos = mapToItem(parent, mouse.x, mouse.y)
                                dragStartY = parentPos.y
                                marginStartY = window.pageMargins.bottom
                            }
                            
                            onPositionChanged: (mouse) => {
                                if (pressed) {
                                    let parentPos = mapToItem(parent, mouse.x, mouse.y)
                                    let deltaY = parentPos.y - dragStartY
                                    let deltaCm = deltaY / parent.cmPixels
                                    let newMargin = marginStartY - deltaCm
                                    
                                    let pHeightCm = window.pageSizeMm(mainToolbar.pageSizeCombo.currentText).h / 10
                                    let maxMargin = pHeightCm - window.pageMargins.top - 1.0
                                    newMargin = Math.max(0.0, Math.min(newMargin, maxMargin))
                                    
                                    let margins = {
                                        left: window.pageMargins.left,
                                        right: window.pageMargins.right,
                                        top: window.pageMargins.top,
                                        bottom: newMargin
                                    }
                                    window.pageMargins = margins
                                }
                            }
                        }
                    }

                    // Plain Item as the RowLayout-managed cell (Layout.fillWidth/
                    // fillHeight live here), so scrollView and the empty-canvas
                    // overlay below can both anchors.fill: parent as ordinary
                    // siblings without QML's "anchors on a layout-managed item is
                    // undefined behavior" warning (confirmed hit when scrollView
                    // itself carried both the Layout.* properties and was the
                    // anchor target for a RowLayout-level sibling).
                    Item {
                        Layout.fillWidth: true
                        Layout.fillHeight: true

                    ScrollView {
                        id: scrollView
                        anchors.fill: parent
                        clip: true
                        ScrollBar.horizontal.policy: ScrollBar.AlwaysOn
                        ScrollBar.vertical.policy: ScrollBar.AlwaysOn
                        // ScrollView infers its scrollable extent from the content item's
                        // IMPLICIT size, but workspaceContent is a plain Item whose explicit
                        // width/height (below) never becomes implicitWidth/implicitHeight
                        // automatically (unlike Column/Row/Layout types). Without this, ScrollView
                        // thought content was always exactly viewport-sized -- the scrollbar thumb
                        // filled the whole track (nothing to grab) and wheel/drag scrolling did
                        // nothing. Same bug class already hit and fixed in components/
                        // StatusPlaceholder.qml; see also centerViewOnPage()'s own comment above,
                        // which worked around half of this same problem for a different call site.
                        contentWidth: workspaceContent.width
                        contentHeight: workspaceContent.height

                        Item {
                            id: workspaceContent
                            width: Math.max(scrollView.availableWidth, (docRect.width * window.zoomLevel) + 100)
                            height: Math.max(scrollView.availableHeight, (docRect.height * window.zoomLevel) + 100)

                            Rectangle {
                                anchors.fill: parent
                                color: Theme.workspaceBackground
                            }

                            Rectangle {
                                id: docRect
                                width: {
                                    switch(mainToolbar.pageSizeCombo.currentText) {
                                        case "A4": return 794;
                                        case "A3": return 1123;
                                        case "A5": return 559;
                                        case "Letter": return 816;
                                        default: return 794;
                                    }
                                }
                                height: {
                                    switch(mainToolbar.pageSizeCombo.currentText) {
                                        case "A4": return 1123;
                                        case "A3": return 1587;
                                        case "A5": return 794;
                                        case "Letter": return 1056;
                                        default: return 1123;
                                    }
                                }
                                color: Theme.surface
                                scale: window.zoomLevel
                                transformOrigin: Item.TopLeft
                                
                                x: Math.max(50, (parent.width - (width * window.zoomLevel)) / 2)
                                y: 50

                                Repeater {
                                    id: canvasRepeater
                                    model: DocumentManager.docIds
                                    anchors.fill: parent

                                    ChemCanvas {
                                        id: canvasInstance
                                        required property int modelData
                                        anchors.fill: parent
                                        visible: modelData === DocumentManager.activeDocId
                                        sketch: DocumentManager.documentFor(modelData)

                                        onIsDirtyChanged: {
                                            var d = window.dirtyDocs
                                            d[modelData] = isDirty
                                            window.dirtyDocs = d
                                        }
                                        Component.onCompleted: {
                                            if (modelData === DocumentManager.activeDocId) window.activeCanvas = canvasInstance
                                        }

                                        onAtomPropertiesRequested: (atomId) => {
                                            if (canvasInstance.sketch) canvasInstance.sketch.requestAtomProperties(atomId)
                                        }

                                        onTextEditRequested: (textId, content, chemX, chemY, isBold, isItalic) => {
                                            taskDialogsGroup.textDialog.textId = textId
                                            taskDialogsGroup.textDialog.chemX = chemX
                                            taskDialogsGroup.textDialog.chemY = chemY
                                            taskDialogsGroup.textDialog.inputText = content
                                            taskDialogsGroup.textDialog.isBold = isBold
                                            taskDialogsGroup.textDialog.isItalic = isItalic
                                            taskDialogsGroup.textDialog.open()
                                        }
                                        onImageInsertRequested: (cx, cy) => {
                                            fileDialogsGroup.imageFileDialog.chemX = cx
                                            fileDialogsGroup.imageFileDialog.chemY = cy
                                            fileDialogsGroup.imageFileDialog.open()
                                        }

                                        Connections {
                                            target: canvasInstance.sketch
                                            function onStructureReady(reqId, data) {
                                                if (reqId === "autosave") {
                                                    // This is for legacy compatibility - new autosave uses requestProjectSerialization
                                                    // which is handled by the structureResponseRoutes above
                                                    const appData = StandardPaths.writableLocation(StandardPaths.AppDataLocation)
                                                    fileIO.write(appData + "/recovery/autosave_" + canvasInstance.modelData + ".qtsketch", data)
                                                } else if (reqId === "autosave_template") {
                                                    window._finishTemplateSave(data)
                                                }
                                            }
                                        }
                                    }
                                }

                                Connections {
                                    target: DocumentManager
                                    function onActiveDocIdChanged() {
                                        for (var i = 0; i < canvasRepeater.count; i++) {
                                            var item = canvasRepeater.itemAt(i)
                                            if (item && item.modelData === DocumentManager.activeDocId) {
                                                window.activeCanvas = item
                                            }
                                        }
                                    }
                                }

                            }
                        }
                    }

                    // A true sibling of scrollView (not a child of docRect, and not
                    // placed inside scrollView's own content either, to avoid
                    // relying on how ScrollView handles extra non-Flickable
                    // children), both anchors.fill: parent within the wrapper Item
                    // above. docRect is the full, zoomable/scrollable page (794x1123
                    // at 100% for A4) -- centering on IT places the hint at the
                    // page's true center, which sits below the visible viewport at
                    // default scroll/zoom (confirmed empirically: a debug marker
                    // centered on docRect rendered well below the on-screen page
                    // area). Filling the wrapper instead overlays it exactly on the
                    // visible viewport, independent of zoom/scroll, and being
                    // declared after scrollView here means it paints on top at the
                    // same default z.
                    Item {
                        anchors.fill: parent
                        StatusPlaceholder {
                            anchors.centerIn: parent
                            visible: window.activeCanvas && propPanel.molAtoms === 0
                            mode: "empty"
                            message: "Empty Canvas"
                            detail: "Select a tool on the left to start drawing,\nor drop a molecule file here"
                        }
                    }
                    } // end wrapper Item
                }
            }

            // Thin separator (hidden when the docked panel is collapsed/floating --
            // there's nothing for it to separate from in that state).
            Rectangle {
                Layout.preferredWidth: 1
                Layout.fillHeight: true
                color: Theme.outline
                visible: uiSettings.panelState === "docked"
            }

            PropertyPanel {
                id: propPanel
                Layout.preferredWidth: uiSettings.panelState === "docked" ? Theme.propertyPanelWidth : 0
                Layout.fillHeight: true
                visible: uiSettings.panelState === "docked"
                canvas: activeCanvas
            }

            // Slim re-expand tab, shown only when the panel isn't docked (covers
            // both "collapsed" and, once Task 4 lands, "floating").
            Rectangle {
                Layout.preferredWidth: 20
                Layout.fillHeight: true
                visible: uiSettings.panelState !== "docked"
                color: Theme.surface

                Text {
                    anchors.centerIn: parent
                    text: "◀"
                    font.pixelSize: Theme.fontSizeCaption
                    color: Theme.textSecondary
                }

                MouseArea {
                    anchors.fill: parent
                    cursorShape: Qt.PointingHandCursor
                    onClicked: uiSettings.panelState = "docked"
                }
            }
        }

    Window {
        id: floatingPropertyWindow
        title: "Properties"
        width: 280
        height: 420
        visible: uiSettings.panelState === "floating"
        flags: Qt.Window

        // Restore the saved position once, before the window first shows --
        // not reactively on panelState change, so a later re-float after a
        // dock/undock cycle doesn't re-run this (the window keeps whatever
        // position the user last dragged it to within the same session).
        Component.onCompleted: {
            if (uiSettings.floatingPanelX < 0 || uiSettings.floatingPanelY < 0) return
            let fits = false
            for (const screen of Qt.application.screens) {
                if (uiSettings.floatingPanelX >= screen.virtualX &&
                    uiSettings.floatingPanelX < screen.virtualX + screen.width &&
                    uiSettings.floatingPanelY >= screen.virtualY &&
                    uiSettings.floatingPanelY < screen.virtualY + screen.height) {
                    fits = true
                    break
                }
            }
            if (fits) {
                x = uiSettings.floatingPanelX
                y = uiSettings.floatingPanelY
            } else {
                // Saved monitor is gone (disconnected, resolution changed,
                // etc.) -- fall back to a position relative to the main
                // window instead of restoring stale, possibly-unreachable
                // coordinates.
                x = window.x + 40
                y = window.y + 40
            }
        }

        onXChanged: if (visible) { uiSettings.floatingPanelX = x }
        onYChanged: if (visible) { uiSettings.floatingPanelY = y }

        // Closing via the OS [X] must not strand panelState as "floating"
        // with the window gone and the docked slot still collapsed to 0 --
        // that would make the panel unreachable. Dock instead of closing.
        onClosing: function(closeEvent) {
            closeEvent.accepted = false
            uiSettings.panelState = "docked"
        }

        PropertyPanel {
            id: floatingPropPanel
            anchors.fill: parent
            canvas: activeCanvas
            isFloating: true
        }
    }

        // Bottom bar: quick-atom picker + status/zoom controls, merged into
        // one row (was two stacked bars costing ~80px combined; see
        // docs/superpowers/specs/2026-10-06-layout-space-maximization-design.md).
        // Tall enough to actually contain a standard-height control plus its
        // padding; at a tighter fixed height the controls inside were
        // squeezed, which is what made the zoom handle read as oversized.
        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: Theme.controlHeight + 8
            color: Theme.surface

            Rectangle {
                anchors.top: parent.top
                width: parent.width
                height: 1
                color: Theme.outline
            }

            RowLayout {
                anchors {
                    left: parent ? parent.left : undefined
                    leftMargin: 10
                    verticalCenter: parent ? parent.verticalCenter : undefined
                }
                spacing: Theme.spacingSmall

                Repeater {
                    model: [
                        { id: "ATOM_C",  icon: "C",  tip: "Carbon" },
                        { id: "ATOM_N",  icon: "N",  tip: "Nitrogen" },
                        { id: "ATOM_O",  icon: "O",  tip: "Oxygen" },
                        { id: "ATOM_S",  icon: "S",  tip: "Sulfur" },
                        { id: "ATOM_F",  icon: "F",  tip: "Fluorine" },
                        { id: "ATOM_Cl", icon: "Cl", tip: "Chlorine" },
                        { id: "ATOM_ANY", icon: "period-table.svg", tip: "Periodic Table" },
                        { id: "CHARGE_PLUS",  icon: "charge-plus.svg", tip: "Charge Plus" },
                        { id: "CHARGE_MINUS", icon: "charge-minus.svg", tip: "Charge Minus" }
                    ]
                    delegate: IconCell {
                        id: atomBarCell
                        required property var modelData
                        iconSource: atomBarCell.modelData.icon.indexOf(".svg") !== -1 ? "icons/" + atomBarCell.modelData.icon : ""
                        glyph: atomBarCell.modelData.icon.indexOf(".svg") === -1 ? atomBarCell.modelData.icon : ""
                        hasGlyphColorOverride: atomBarCell.modelData.icon.indexOf(".svg") === -1
                        glyphColor: atomBarCell.modelData.icon === "C" ? Theme.textPrimary : Theme.getElementColor(atomBarCell.modelData.icon)
                        tip: atomBarCell.modelData.tip
                        selected: !!(activeCanvas && activeCanvas.currentTool === atomBarCell.modelData.id)
                        onClicked: {
                            if (!activeCanvas) return
                            if (atomBarCell.modelData.id === "ATOM_ANY") {
                                periodicTablePopup.open()
                            } else if (activeCanvas.currentTool === atomBarCell.modelData.id) {
                                activeCanvas.currentTool = "SELECT"
                            } else {
                                activeCanvas.currentTool = atomBarCell.modelData.id
                            }
                        }
                    }
                }

                Rectangle { width: 1; height: 20; color: Theme.outline; opacity: 0.6 }

                Text {
                    text: activeCanvas ? "Tool: " + ToolLabels.label(activeCanvas.currentTool) : "No document"
                    color: Theme.textSecondary
                    font { pixelSize: Theme.fontSizeBody; family: Theme.fontMono }
                }
                // text: bindings evaluate even while visible is false, so every
                // dereference below stays guarded down to the leaf array.
                Text {
                    visible: !!(activeSketch && activeSketch.primitives && activeSketch.primitives.atoms)
                    text: "Atoms: " + ((activeSketch && activeSketch.primitives && activeSketch.primitives.atoms) ? activeSketch.primitives.atoms.length : 0)
                    color: Theme.textSecondary
                    font { pixelSize: Theme.fontSizeBody; family: Theme.fontMono }
                }
                Text {
                    visible: !!(activeSketch && activeSketch.primitives && activeSketch.primitives.bonds)
                    text: "Bonds: " + ((activeSketch && activeSketch.primitives && activeSketch.primitives.bonds) ? activeSketch.primitives.bonds.length : 0)
                    color: Theme.textSecondary
                    font { pixelSize: Theme.fontSizeBody; family: Theme.fontMono }
                }
                Text {
                    visible: !!(activeSketch && activeSketch.selection && activeSketch.selection.atom_ids)
                    // Split atoms/bonds instead of one summed count -- the summed
                    // form ("Selected: 29") contradicted the Selection Summary
                    // panel ("14 atoms, 15 bonds selected") for the same state.
                    text: {
                        const sel = activeSketch ? activeSketch.selection : null
                        if (!sel) return ""
                        const na = sel.atom_ids ? sel.atom_ids.length : 0
                        const nb = sel.bond_ids ? sel.bond_ids.length : 0
                        const rest = Selection.totalCount(activeSketch) - na - nb
                        let parts = "Selected: " + na + "a " + nb + "b"
                        if (rest > 0) parts += " +" + rest
                        return parts
                    }
                    color: Theme.textSecondary
                    font { pixelSize: Theme.fontSizeBody; family: Theme.fontMono }
                }

                Rectangle {
                    property int issueCount: {
                        if (!activeSketch || !activeSketch.primitives) return 0;
                        let count = 0;
                        if (activeSketch.primitives.atoms) {
                            for (let i = 0; i < activeSketch.primitives.atoms.length; i++) {
                                if (activeSketch.primitives.atoms[i].checkWarning) count++;
                            }
                        }
                        if (activeSketch.primitives.bonds) {
                            for (let j = 0; j < activeSketch.primitives.bonds.length; j++) {
                                if (activeSketch.primitives.bonds[j].checkWarning) count++;
                            }
                        }
                        return count;
                    }
                    visible: issueCount > 0
                    width: issueRow.implicitWidth + 12
                    height: 20
                    radius: 10
                    color: Theme.error

                    RowLayout {
                        id: issueRow
                        anchors.centerIn: parent
                        spacing: 4
                        Text {
                            // Segoe Fluent Icons: Warning (U+E7BA). Built from a
                            // char code rather than pasted inline — private-use
                            // codepoints do not survive every editing path.
                            text: String.fromCharCode(0xE7BA)
                            color: Theme.badgeText
                            font { family: Theme.fontIcons; pixelSize: 12 }
                        }
                        Text {
                            text: parent.parent.issueCount + (parent.parent.issueCount === 1 ? " issue" : " issues")
                            color: Theme.badgeText
                            font { pixelSize: Theme.fontSizeCaption; family: Theme.fontFamily }
                        }
                    }
                    MouseArea {
                        anchors.fill: parent
                        cursorShape: Qt.PointingHandCursor
                        onClicked: {
                            // Re-expand the Properties panel if it's collapsed/floating --
                            // wired in Task 3/4; a no-op here if already docked.
                            if (typeof uiSettings !== "undefined" && uiSettings.panelState !== "docked") {
                                uiSettings.panelState = "docked"
                            }
                        }
                    }
                }
            }

            RowLayout {
                anchors {
                    right: parent ? parent.right : undefined
                    rightMargin: 16
                    verticalCenter: parent ? parent.verticalCenter : undefined
                }
                spacing: 8

                CheckBox {
                    text: "Show Hydrogens"
                    checked: activeCanvas ? activeCanvas.showExplicitH : false
                    onCheckedChanged: if (activeCanvas) activeCanvas.showExplicitH = checked
                }

                ToolButton {
                    // ViewAll (large cells) when compact is active, GridView
                    // (dense cells) when comfortable is — the icon shows the
                    // layout the click switches TO, matching the tooltip.
                    text: Theme.compactDensity ? "" : ""
                    font { family: Theme.fontIcons; pixelSize: Theme.iconSize }
                    implicitWidth: Theme.controlHeight
                    implicitHeight: Theme.controlHeight
                    Accessible.name: Theme.compactDensity ? "Switch to comfortable density" : "Switch to compact density"
                    onClicked: {
                        Theme.compactDensity = !Theme.compactDensity
                        uiSettings.compactDensity = Theme.compactDensity
                    }
                    ToolTip.text: Theme.compactDensity ? "Switch to comfortable density" : "Switch to compact density"
                    ToolTip.visible: hovered
                    ToolTip.delay: 500
                }

                Text {
                    text: Math.round(window.zoomLevel * 100) + "%"
                    color: Theme.textSecondary
                    font { pixelSize: Theme.fontSizeBody; bold: true; family: Theme.fontMono }
                    MouseArea {
                        anchors.fill: parent
                        cursorShape: Qt.PointingHandCursor
                        onClicked: window.zoomLevel = 1.0
                    }
                }

                ToolButton {
                    // Segoe Fluent Icons ZoomOut (U+E71F). The SVG version drew
                    // nothing under the Fluent style even though it reserved its
                    // space, so the glyph font is used here as it is for the
                    // density toggle beside it.
                    text: String.fromCharCode(0xE71F)
                    font { family: Theme.fontIcons; pixelSize: Theme.iconSize }
                    implicitWidth: Theme.controlHeight
                    implicitHeight: Theme.controlHeight
                    Accessible.name: "Zoom out"
                    onClicked: window.zoomLevel = Math.max(0.1, window.zoomLevel - 0.1)
                }

                Slider {
                    id: zoomSlider
                    from: 0.1
                    to: 3.0
                    Layout.preferredWidth: 120
                    Layout.preferredHeight: Theme.controlHeight
                    onMoved: window.zoomLevel = value
                    Binding on value {
                        value: window.zoomLevel
                        restoreMode: Binding.RestoreBinding
                    }
                }

                ToolButton {
                    // Segoe Fluent Icons ZoomIn (U+E8A3)
                    text: String.fromCharCode(0xE8A3)
                    font { family: Theme.fontIcons; pixelSize: Theme.iconSize }
                    implicitWidth: Theme.controlHeight
                    implicitHeight: Theme.controlHeight
                    Accessible.name: "Zoom in"
                    onClicked: window.zoomLevel = Math.min(3.0, window.zoomLevel + 0.1)
                }
            }
        }
    }

    property url pendingSdfBatchUrl: ""
    property string _pendingBatchAction: "scaffold"
    property var _pendingBatchLabels: []
    property var _pendingBatchGridMolfiles: []
    property string pendingSimilarityRefMolfile: ""
    property var _rgroupPerMoleculeResults: []

    FileIO {
        id: fileIO
    }

    SdfRecordPicker {
        id: sdfRecordPicker
        onRecordChosen: function(index) {
            if (activeSketch) activeSketch.sendCommand("loadSdfBatchRecord", [index])
        }
        onFindScaffoldRequested: {
            window.isProcessing = true
            _pendingBatchAction = "scaffold"
            if (activeSketch) activeSketch.sendCommand("getSdfBatchMolfiles")
            sdfRecordPicker.close()
        }
        onDecomposeRequested: {
            window.isProcessing = true
            _pendingBatchAction = "decompose"
            if (activeSketch) activeSketch.sendCommand("getSdfBatchMolfiles")
            sdfRecordPicker.close()
        }
        onDecomposePerMoleculeRequested: {
            window.isProcessing = true
            _pendingBatchAction = "decompose_permolecule"
            if (activeSketch) activeSketch.sendCommand("getSdfBatchMolfiles")
            sdfRecordPicker.close()
        }
        onRankBySimilarityRequested: {
            window.isProcessing = true
            _pendingBatchAction = "similarity"
            window.pendingSimilarityRefMolfile = activeCanvas ? activeCanvas.getMolfile() : ""
            if (activeSketch) activeSketch.sendCommand("getSdfBatchMolfiles")
            sdfRecordPicker.close()
        }
        onAlignToScaffoldRequested: {
            window.isProcessing = true
            _pendingBatchAction = "align"
            if (activeSketch) activeSketch.sendCommand("getSdfBatchMolfiles")
        }
        onExportGridRequested: {
            _pendingBatchAction = "export_grid"
            if (activeSketch) activeSketch.sendCommand("getSdfBatchMolfiles")
            sdfRecordPicker.close()
        }
        onExportBatchFileRequested: {
            _pendingBatchAction = "export_file"
            if (activeSketch) activeSketch.sendCommand("getSdfBatchMolfiles")
        }
    }

    SubstructureSearchPopup {
        id: smartsSearchPopup
        onSearchRequested: function(smarts) {
            if (!activeSketch) return
            window.pendingSmartsQuery = smarts
            activeSketch.requestSerialize("smarts_search")
        }
        onClearRequested: {
            if (activeSketch) activeSketch.sendCommand("selectSubstructureMatches", ["{\"matches\":[]}"])
        }
    }
    property string pendingSmartsQuery: ""
    property string pendingSimilarityRef: ""
    property double pendingIonizePh: 7.4

    function pageSizeMm(name) {
        switch (name) {
            case "A4": return { w: 210, h: 297 }
            case "A3": return { w: 297, h: 420 }
            case "A5": return { w: 148, h: 210 }
            case "Letter": return { w: 215.9, h: 279.4 }
            default: return { w: 210, h: 297 }
        }
    }

    function printToPdf() {
        if (!activeCanvas) return
        fileDialogsGroup.pdfSaveDialog.open()
    }

    function loadFromFile(fileUrl) {
        const fileStrEarly = fileUrl.toString().toLowerCase()
        if (fileStrEarly.endsWith(".rdf")) {
            pendingSdfBatchUrl = fileUrl
            window.isProcessing = true
            indigoSvc.parseRdfBatch(fileUrl)
            return
        } else if (fileStrEarly.endsWith(".smi") || fileStrEarly.endsWith(".smiles")) {
            pendingSdfBatchUrl = fileUrl
            window.isProcessing = true
            indigoSvc.parseIndigoBatchFile(fileUrl, "smiles")
            return
        } else if (fileStrEarly.endsWith(".cml")) {
            pendingSdfBatchUrl = fileUrl
            window.isProcessing = true
            indigoSvc.parseIndigoBatchFile(fileUrl, "cml")
            return
        } else if (fileStrEarly.endsWith(".cdx")) {
            pendingSdfBatchUrl = fileUrl
            window.isProcessing = true
            indigoSvc.parseIndigoBatchFile(fileUrl, "cdx")
            return
        } else if (fileStrEarly.endsWith(".rxn")) {
            const data = fileIO.read(fileUrl)
            if (data !== "" && activeSketch) {
                if (activeCanvas) activeCanvas.clearCanvas()
                const ok = activeSketch.importReaction(data)
                if (ok) {
                    window.centerViewOnPage()
                    setDocFile(DocumentManager.activeDocId, fileUrl)
                } else {
                    messageDialogsGroup.workerErrorDialog.errorText = "Open Reaction: the file does not contain a valid reaction. The canvas was cleared before the import attempt -- press Ctrl+Z to restore your previous content."
                    messageDialogsGroup.workerErrorDialog.severe = false
                    messageDialogsGroup.workerErrorDialog.open()
                }
            }
            return
        }
        const data = fileIO.read(fileUrl)
        if (data === "") return
        const fileStr = fileUrl.toString().toLowerCase()
        if (fileStr.endsWith(".fasta") || fileStr.endsWith(".fa")) {
            biopolymerDialog.open()
            window.isProcessing = true
            indigoSvc.loadBioFasta(data, "DNA")
        } else if (fileStr.endsWith(".helm")) {
            window.isProcessing = true
            indigoSvc.loadBioHelm(data)
        } else if (fileStr.endsWith(".idt")) {
            window.isProcessing = true
            indigoSvc.loadBioIdt(data)
        } else {
            let fmt = "mol"
            if (fileStr.endsWith(".sdf")) fmt = "sdf"
            else if (fileStr.endsWith(".ket")) fmt = "ket"
            else if (fileStr.endsWith(".qtsketch")) fmt = "qtsketch"
            
            if (fmt === "sdf" && activeSketch) {
                pendingSdfBatchUrl = fileUrl
                activeSketch.sendCommand("deserializeSdfBatch", [data])
            } else if (fmt === "qtsketch" && activeSketch) {
                // Use the new project load method - pass the file path as string
                const filePath = fileUrl.toString()
                activeSketch.loadProject(filePath)
                window.centerViewOnPage()
                setDocFile(DocumentManager.activeDocId, fileUrl)
            } else if (activeCanvas) {
                activeCanvas.loadStructure(fmt, data)
                window.centerViewOnPage()
                setDocFile(DocumentManager.activeDocId, fileUrl)
            }
        }
    }

    DropArea {
        anchors.fill: parent
        enabled: true
        onDropped: (drop) => {
            if (drop.hasUrls && drop.urls.length > 0) {
                loadFromFile(drop.urls[0])
                drop.accept()
            }
        }
    }

    PeriodicTablePopup {
        id: periodicTablePopup
        onElementSelected: (label) => {
            activeCanvas.currentTool = "ATOM_" + label
        }
    }

    Shortcut {
        sequence: "Ctrl+I"
        onActivated: taskDialogsGroup.smilesDialog.open()
    }

    Shortcut {
        sequence: "Ctrl+B"
        onActivated: biopolymerDialog.open()
    }

    Shortcut {
        sequence: "Ctrl+R"
        onActivated: rgroupPanel.open()
    }

    Shortcut {
        sequence: "Ctrl+N"
        onActivated: DocumentManager.addDocument()
    }

    Shortcut {
        sequence: "Ctrl+S"
        onActivated: window.saveActive(false)
    }

    Shortcut {
        sequence: "Ctrl+Shift+S"
        onActivated: window.saveActive(true)  // Save As: always prompt
    }

    Shortcut {
        sequence: "Ctrl+O"
        onActivated: fileDialogsGroup.openDialog.open()
    }

    Shortcut {
        sequence: "Ctrl+Shift+Z"
        onActivated: if (activeCanvas) activeCanvas.redo()
    }

    Shortcut {
        sequence: "Ctrl+W"
        onActivated: {
            if (DocumentManager.docIds.length > 1) {
                if (window.dirtyDocs[DocumentManager.activeDocId]) {
                    window._pendingCloseDocId = DocumentManager.activeDocId
                    messageDialogsGroup.unsavedChangesDialog.open()
                } else {
                    DocumentManager.closeDocument(DocumentManager.activeDocId)
                }
            }
        }
    }

    RGroupPanel {
        id: rgroupPanel
        sketch: window.activeSketch
        canvas: window.activeCanvas
    }

    BiopolymerSequenceView {
        id: biopolymerSeqView
        sketch: window.activeSketch
    }

    BiopolymerDialog {
        id: biopolymerDialog

        onPreviewSequenceViewRequested: function(seqType, text) {
            if (activeSketch) activeSketch.sendCommand("bioBuildSequenceView", [text, seqType])
        }

        onLoadRequested: function(format, seqType, text) {
            window.isProcessing = true
            if      (format === "Sequence") indigoSvc.loadBioSequence(text, seqType)
            else if (format === "FASTA")    indigoSvc.loadBioFasta(text, seqType)
            else if (format === "HELM")     indigoSvc.loadBioHelm(text)
            else if (format === "IDT")      indigoSvc.loadBioIdt(text)
            else if (format === "AxoLabs")  indigoSvc.loadBioAxoLabs(text)
        }

        onExportRequested: function(format) {
            if      (format === "Sequence") activeSketch.requestSerialize("bio_seq")
            else if (format === "FASTA")    activeSketch.requestSerialize("bio_fasta")
            else if (format === "HELM")     activeSketch.requestSerialize("bio_helm")
            else if (format === "IDT")      activeSketch.requestSerialize("bio_idt")
            else if (format === "AxoLabs")  activeSketch.requestSerialize("bio_axolabs")
        }
    }

    AtomPropertiesDialog {
        id: atomPropsDialog
        onPropertiesChanged: function(atomId, isotope, radical, valence) {
            if (activeCanvas) activeCanvas.applyPropertyChange("atomProps", atomId, { isotope: isotope, radical: radical, valence: valence })
        }
    }

    FileDialogs {
        id: fileDialogsGroup
        win: window
    }
    TaskDialogs {
        id: taskDialogsGroup
        win: window
    }

    MessageDialogs {
        id: messageDialogsGroup
        win: window
    }

    // TaskDialog rather than a raw Dialog: TaskDialog already supplies
    // modal:true + anchors.centerIn:parent (components/TaskDialog.qml), which
    // this dialog was re-deriving by hand via explicit x/y bindings on
    // window.width/height -- every other dialog in the app goes through the
    // shared base, and this was the one exception.
    TaskDialog {
        id: userTemplateDialog
        title: "Save Template"
        standardButtons: Dialog.Ok | Dialog.Cancel
        width: Theme.dialogWidthSmall

        ColumnLayout {
            width: parent.width
            Label { text: "Template Name:" }
            TextField {
                id: templateNameInput
                Layout.fillWidth: true
                focus: true
                onAccepted: userTemplateDialog.accept()
            }
        }
        
        onAccepted: {
            if (window.activeSketch) {
                window._pendingTemplateName = templateNameInput.text
                window.activeSketch.requestSelectionStructure("autosave_template")
            }
        }
        onOpened: templateNameInput.text = ""
    }

    function saveSelectionAsTemplate() {
        if (!activeSketch) return
        userTemplateDialog.open()
    }

    property string _pendingTemplateName: ""

    function _finishTemplateSave(data) {
        if (!data || data === "") return
        var name = _pendingTemplateName
        if (!name) name = "Unnamed"
        window._userTemplates.push({ name: name, data: data })
        window.saveUserTemplates()
        // force property update
        var tmp = window._userTemplates
        window._userTemplates = []
        window._userTemplates = tmp
    }

    Popup {
        id: imageConfirmBanner
        x: (window.width - width) / 2
        y: window.height - height - 30
        width: 400
        height: 50
        modal: false
        focus: false
        closePolicy: Popup.NoAutoClose

        background: Rectangle {
            color: Theme.surface
            radius: 4
            border.color: Theme.outline
            border.width: 1
        }

        function showBanner(warnings) {
            bannerText.text = "Recognized structure inserted (" + warnings + " warning" + (warnings === 1 ? "" : "s") + ")."
            bannerTimer.restart()
            open()
        }

        RowLayout {
            anchors.fill: parent
            anchors.margins: 8
            spacing: 8
            Text {
                id: bannerText
                color: Theme.textPrimary
                font.pixelSize: Theme.fontSizeBody
                Layout.fillWidth: true
            }
            Button {
                text: "Keep"
                onClicked: { imageConfirmBanner.close(); bannerTimer.stop() }
            }
            Button {
                text: "Undo & use image instead"
                onClicked: {
                    imageConfirmBanner.close()
                    bannerTimer.stop()
                    if (window.activeCanvas) window.activeCanvas.undo()
                    const dataUri = fileIO.readImageAsDataUri(window.fileDialogsGroup.imageFileDialog.pendingFileUrl)
                    if (dataUri !== "" && window.activeSketch) {
                        window.activeSketch.addImage(dataUri, window.fileDialogsGroup.imageFileDialog.chemX, window.fileDialogsGroup.imageFileDialog.chemY, 1.5, 1.5)
                    }
                }
            }
        }

        Timer {
            id: bannerTimer
            interval: 6000
            repeat: false
            onTriggered: imageConfirmBanner.close()
        }
    }
}
