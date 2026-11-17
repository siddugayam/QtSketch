import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtCore

Rectangle {
    id: root

    property string currentTool: "SELECT"
    property var win: undefined
    // Bound by MainWindow to the currently-active document's V8Process instance.
    property var sketch: null

    signal toolSelected(string toolId)

    color: Theme.surface

    // ── Collapsible section state (persisted across restarts) ───────────────
    property bool editOpen: true
    property bool bondsOpen: true
    property bool ringsOpen: true
    property bool groupsOpen: true
    property bool templatesOpen: true
    property bool userTemplatesOpen: false

    Settings {
        category: "toolPanel"
        property alias editOpen: root.editOpen
        property alias bondsOpen: root.bondsOpen
        property alias ringsOpen: root.ringsOpen
        property alias groupsOpen: root.groupsOpen
        property alias templatesOpen: root.templatesOpen
    }

    // Auto-expand the section owning a tool activated by other means (single-key
    // shortcut, periodic-table pick). Only ever opens sections — never collapses.
    onCurrentToolChanged: {
        if (currentTool === "SELECT" || currentTool === "SELECT_FRAGMENT" || currentTool === "SELECT_LASSO" ||
            currentTool === "HAND" || currentTool === "ERASE" || currentTool === "TEXT" || currentTool === "TEXT_NUMBER")
            editOpen = true
        else if (currentTool.indexOf("BOND_") === 0 || currentTool === "CHAIN")
            bondsOpen = true
        else if (currentTool.indexOf("TEMPLATE_") === 0 || currentTool === "LIB_Pyridine")
            ringsOpen = true
        else if (currentTool.indexOf("FG_") === 0 || currentTool.indexOf("SS_") === 0)
            groupsOpen = true
        else if (currentTool.indexOf("LIB_") === 0)
            templatesOpen = true
    }

    // ── Reusable pieces ──────────────────────────────────────────────────────
    component SectionHeader: Rectangle {
        id: header
        property string title
        property bool expanded: true
        signal toggled()
        Layout.fillWidth: true
        implicitHeight: 22
        radius: 3
        color: headerMouse.containsMouse ? Theme.hover : "transparent"
        Accessible.role: Accessible.Button
        Accessible.name: header.title
        Behavior on color { ColorAnimation { duration: 150 } }

        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: 4
            spacing: 4
            Text {
                // Segoe Fluent Icons ChevronDown / ChevronRight, written as
                // \u escapes rather than literal private-use characters so the
                // source stays ASCII and the glyphs cannot be lost in transit.
                text: header.expanded ? Theme.icon(0xE70D, 0xF0140) : Theme.icon(0xE76C, 0xF0142)
                color: Theme.textSecondary
                font { family: Theme.fontIcons; pixelSize: 10 }
            }
            Text {
                Layout.fillWidth: true
                text: header.title
                color: Theme.textSecondary
                elide: Text.ElideRight
                // Fluent section headers are sentence-case Body Strong, not
                // letterspaced small caps. Dropping the tracking and the caps
                // also buys back the width that was truncating "Templates".
                font {
                    pixelSize: Theme.fontSizeLabel
                    weight: Font.DemiBold
                    family: Theme.fontFamily
                }
            }
        }
        MouseArea {
            id: headerMouse
            anchors.fill: parent
            hoverEnabled: true
            cursorShape: Qt.PointingHandCursor
            onClicked: header.toggled()
        }
    }

    // Scrollable so every tool stays reachable when the window is shorter than
    // the stacked tool sections (~700px of content).
    ScrollView {
        id: toolScroll
        anchors.fill: parent
        padding: Theme.marginSmall
        clip: true
        contentWidth: availableWidth
        ScrollBar.horizontal.policy: ScrollBar.AlwaysOff

        ColumnLayout {
        width: toolScroll.availableWidth
        spacing: Theme.spacingSmall

        // ── Edit section ─────────────────────────────────────────────────────
        SectionHeader {
            title: "Edit"
            expanded: root.editOpen
            onToggled: root.editOpen = !root.editOpen
        }
        Item {
            Layout.fillWidth: true
            Layout.preferredHeight: root.editOpen ? editGrid.implicitHeight : 0
            Behavior on Layout.preferredHeight { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }
            clip: true
            opacity: root.editOpen ? 1 : 0
            Behavior on opacity { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }

            GridLayout {
                id: editGrid
                anchors.top: parent.top
                anchors.horizontalCenter: parent.horizontalCenter
                columns: 3
            columnSpacing: Theme.toolGridGap
            rowSpacing: Theme.toolGridGap
            Repeater {
                model: [
                    { id: "SELECT", icon: "select.svg", tip: "Selection tool (S)" },
                    { id: "SELECT_FRAGMENT", icon: "select-fragment.svg", tip: "Fragment Selection tool" },
                    { id: "SELECT_RING", icon: "Ring", tip: "Ring Selection tool" },
                    { id: "SELECT_CHAIN", icon: "Chain", tip: "Chain Selection tool" },
                    { id: "SELECT_LASSO", icon: "select-lasso.svg", tip: "Lasso select tool - freeform selection" },
                    { id: "HAND", icon: "hand.svg", tip: "Hand tool - drag to pan canvas (H)" },
                    { id: "ERASE", icon: "erase.svg", tip: "Erase tool (E)" },
                    { id: "TEXT", icon: "text.svg", tip: "Text annotation (click canvas)" },
                    { id: "TEXT_NUMBER", icon: "1a", tip: "Auto-numbering annotation" },
                    { id: "IMAGE", icon: "add-image.svg", tip: "Insert image" }
                ]
                delegate: IconCell {
                    required property var modelData
                    Layout.alignment: Qt.AlignHCenter
                    iconSource: modelData.icon.indexOf(".svg") !== -1 ? "icons/" + modelData.icon : ""
                    glyph: modelData.icon.indexOf(".svg") === -1 ? modelData.icon : ""
                    tip: modelData.tip
                    selected: root.currentTool === modelData.id
                    onClicked: root.toolSelected(modelData.id)
                }
            }
        }
        }

        // ── Bonds section ────────────────────────────────────────────────────
        SectionHeader {
            title: "Bonds"
            expanded: root.bondsOpen
            onToggled: root.bondsOpen = !root.bondsOpen
        }
        Item {
            Layout.fillWidth: true
            Layout.preferredHeight: root.bondsOpen ? bondsGrid.implicitHeight : 0
            Behavior on Layout.preferredHeight { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }
            clip: true
            opacity: root.bondsOpen ? 1 : 0
            Behavior on opacity { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }

            GridLayout {
                id: bondsGrid
                anchors.top: parent.top
                anchors.horizontalCenter: parent.horizontalCenter
                columns: 3
            columnSpacing: Theme.toolGridGap
            rowSpacing: Theme.toolGridGap
            Repeater {
                model: [
                    { id: "BOND_1", icon: "single_bond.svg", tip: "Single bond (B)" },
                    { id: "BOND_2", icon: "double_bond.svg", tip: "Double bond" },
                    { id: "BOND_3", icon: "triple_bond.svg", tip: "Triple bond" },
                    { id: "BOND_UP", icon: "up_bond.svg", tip: "Wedge bond" },
                    { id: "BOND_DOWN", icon: "down_bond.svg", tip: "Dash bond" },
                    { id: "BOND_UPDOWN", icon: "updown_bond.svg", tip: "Wavy bond" },
                    { id: "CHAIN", icon: "chain.svg", tip: "Chain (drag to lay a carbon chain)" }
                ]
                delegate: IconCell {
                    required property var modelData
                    Layout.alignment: Qt.AlignHCenter
                    iconSource: modelData.icon.indexOf(".svg") !== -1 ? "icons/" + modelData.icon : ""
                    glyph: modelData.icon.indexOf(".svg") === -1 ? modelData.icon : ""
                    tip: modelData.tip
                    selected: root.currentTool === modelData.id
                    onClicked: root.toolSelected(modelData.id)
                }
            }
        }
        }

        // ── Rings section ────────────────────────────────────────────────────
        SectionHeader {
            title: "Rings"
            expanded: root.ringsOpen
            onToggled: root.ringsOpen = !root.ringsOpen
        }
        Item {
            Layout.fillWidth: true
            Layout.preferredHeight: root.ringsOpen ? ringsGrid.implicitHeight : 0
            Behavior on Layout.preferredHeight { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }
            clip: true
            opacity: root.ringsOpen ? 1 : 0
            Behavior on opacity { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }

            GridLayout {
                id: ringsGrid
                anchors.top: parent.top
                anchors.horizontalCenter: parent.horizontalCenter
                columns: 3
            columnSpacing: Theme.toolGridGap
            rowSpacing: Theme.toolGridGap
            Repeater {
                model: [
                    // template-N files are in Ketcher's palette order, not ring size:
                    // 0=benzene, 2=hexagon, 3=pentagon, 5=square, 4=triangle
                    { id: "TEMPLATE_BENZENE", icon: "template-0.svg", tip: "Benzene (R)" },
                    { id: "LIB_Pyridine", icon: "Pyr", tip: "Pyridine" },
                    { id: "TEMPLATE_6", icon: "template-2.svg", tip: "Cyclohexane" },
                    { id: "TEMPLATE_5", icon: "template-3.svg", tip: "Cyclopentane" },
                    { id: "TEMPLATE_4", icon: "template-5.svg", tip: "Cyclobutane" },
                    { id: "TEMPLATE_3", icon: "template-4.svg", tip: "Cyclopropane" },
                    { id: "TEMPLATE_7", icon: "template-6.svg", tip: "Cycloheptane" },
                    { id: "TEMPLATE_8", icon: "template-7.svg", tip: "Cyclooctane" }
                ]
                delegate: IconCell {
                    required property var modelData
                    Layout.alignment: Qt.AlignHCenter
                    iconSource: modelData.icon.indexOf(".svg") !== -1 ? "icons/" + modelData.icon : ""
                    glyph: modelData.icon.indexOf(".svg") === -1 ? modelData.icon : ""
                    tip: modelData.tip
                    selected: root.currentTool === modelData.id
                    onClicked: root.toolSelected(modelData.id)
                }
            }
        }
        }

        // ── Groups & Salts section ───────────────────────────────────────────
        SectionHeader {
            title: "Groups"
            expanded: root.groupsOpen
            onToggled: root.groupsOpen = !root.groupsOpen
        }
        Item {
            Layout.fillWidth: true
            Layout.preferredHeight: root.groupsOpen ? groupsGrid.implicitHeight : 0
            Behavior on Layout.preferredHeight { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }
            clip: true
            opacity: root.groupsOpen ? 1 : 0
            Behavior on opacity { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }

            GridLayout {
                id: groupsGrid
                anchors.top: parent.top
                anchors.horizontalCenter: parent.horizontalCenter
                columns: 3
                columnSpacing: Theme.toolGridGap
                rowSpacing: Theme.toolGridGap

                // Functional Groups Tool
                IconCell {
                    Layout.alignment: Qt.AlignHCenter
                    iconSource: "icons/generic-groups.svg"
                    tip: "Functional Groups"
                    selected: root.currentTool.startsWith("FG_")
                    onClicked: {
                        fgList.model = []
                        fgPopup.open()
                        sketch.requestFunctionalGroupsList()
                    }
                }

                // Salts & Solvents Tool
                IconCell {
                    Layout.alignment: Qt.AlignHCenter
                    glyph: "S&S"
                    tip: "Salts & Solvents"
                    selected: root.currentTool.startsWith("SS_")
                    onClicked: {
                        saltsList.model = []
                        ssPopup.open()
                        sketch.requestSaltsAndSolventsList()
                    }
                }
            }
        }

        // ── Templates section (Library launchers) ───────────────
        SectionHeader {
            title: "Templates"
            expanded: root.templatesOpen
            onToggled: root.templatesOpen = !root.templatesOpen
        }
        Item {
            Layout.fillWidth: true
            Layout.preferredHeight: root.templatesOpen ? templatesGrid.implicitHeight : 0
            Behavior on Layout.preferredHeight { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }
            clip: true
            opacity: root.templatesOpen ? 1 : 0
            Behavior on opacity { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }

            GridLayout {
                id: templatesGrid
                anchors.top: parent.top
                anchors.horizontalCenter: parent.horizontalCenter
                columns: 3
            columnSpacing: Theme.toolGridGap
            rowSpacing: Theme.toolGridGap

        // Template Library Tool
        IconCell {
            Layout.alignment: Qt.AlignHCenter
            iconSource: "icons/template-lib.svg"
            tip: "Template Library"
            selected: root.currentTool.startsWith("LIB_")
            onClicked: {
                root._libraryAllItems = []
                libGroupCombo.currentIndex = 0
                libSearchField.text = ""
                libItemList.model = []
                libPopup.open()
                sketch.requestTemplateLibraryList()
            }
        }
        }
        }

        // ── User Templates section ───────────────
        SectionHeader {
            title: "Saved"
            expanded: root.userTemplatesOpen
            onToggled: root.userTemplatesOpen = !root.userTemplatesOpen
        }
        Item {
            Layout.fillWidth: true
            Layout.preferredHeight: root.userTemplatesOpen ? userTemplatesGrid.implicitHeight : 0
            Behavior on Layout.preferredHeight { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }
            clip: true
            opacity: root.userTemplatesOpen ? 1 : 0
            Behavior on opacity { NumberAnimation { duration: Theme.durationMedium; easing.type: Easing.BezierSpline; easing.bezierCurve: Theme.easingDecelerate } }

            GridLayout {
                id: userTemplatesGrid
                anchors.top: parent.top
                anchors.horizontalCenter: parent.horizontalCenter
                columns: 3
                columnSpacing: Theme.toolGridGap
                rowSpacing: Theme.toolGridGap

                Repeater {
                    model: root.win ? root.win._userTemplates : []
                    IconCell {
                        Layout.alignment: Qt.AlignHCenter
                        glyph: modelData.name.substring(0, 2).toUpperCase()
                        tip: modelData.name
                        selected: root.currentTool === "USER_" + index
                        onClicked: {
                            root.toolSelected("USER_" + index)
                        }
                        
                        Rectangle {
                            width: 14
                            height: 14
                            color: Theme.error
                            radius: 7
                            anchors.top: parent.top
                            anchors.right: parent.right
                            // Gated on the enclosing IconCell's own `hovered` (AbstractButton
                            // property), not this Rectangle's own MouseArea.containsMouse --
                            // a hidden Item's MouseArea never receives hover, so gating on its
                            // own containsMouse made the badge permanently unreachable.
                            visible: hovered
                            Text {
                                // Segoe Fluent Icons Dismiss (U+E8BB)
                                text: Theme.icon(0xE8BB, 0xF0156)
                                color: "white"
                                font { family: Theme.fontIcons; pixelSize: 9 }
                                anchors.centerIn: parent
                            }
                            MouseArea {
                                id: mouseAreaForDel
                                anchors.fill: parent
                                hoverEnabled: true
                                onClicked: {
                                    if (root.win) {
                                        root.win._userTemplates.splice(index, 1)
                                        root.win.saveUserTemplates()
                                        var tmp = root.win._userTemplates
                                        root.win._userTemplates = []
                                        root.win._userTemplates = tmp
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        Item { Layout.fillHeight: true }  // Spacer
        }
    }

    property var _saltsListData: []
    property var _libraryAllItems: []
    property var _libraryGroups: []
    property var _thumbnails: ({})

    Connections {
        target: root.sketch
        function onStructureReady(reqId, data) {
            if (reqId === "salts") {
                try {
                    var items = JSON.parse(data)
                    var labels = []
                    for (var i = 0; i < items.length; ++i) labels.push(items[i].label)
                    saltsList.model = labels
                    root._thumbnails = ({})
                    for (var s = 0; s < labels.length; ++s) {
                        sketch.requestTemplateThumbnail(labels[s], "thumb_ss_" + labels[s])
                    }
                } catch (e) {}
            } else if (reqId === "fg_list") {
                try {
                    var fgItems = JSON.parse(data)
                    var fgLabels = []
                    for (var j = 0; j < fgItems.length; ++j) fgLabels.push(fgItems[j].label)
                    fgList.model = fgLabels
                    // Request thumbnails for each FG
                    root._thumbnails = ({})
                    for (var t = 0; t < fgLabels.length; ++t) {
                        sketch.requestTemplateThumbnail(fgLabels[t], "thumb_fg_" + fgLabels[t])
                    }
                } catch (e) {}
            } else if (reqId.indexOf("thumb_fg_") === 0) {
                try {
                    var name = reqId.substring(9)
                    var thumb = JSON.parse(data)
                    var updated = Object.assign({}, root._thumbnails)
                    updated[name] = thumb
                    root._thumbnails = updated
                } catch (e) {}
            } else if (reqId.indexOf("thumb_ss_") === 0) {
                try {
                    var name = reqId.substring(9)
                    var thumb = JSON.parse(data)
                    var updated = Object.assign({}, root._thumbnails)
                    updated[name] = thumb
                    root._thumbnails = updated
                } catch (e) {}
            } else if (reqId.indexOf("thumb_lib_") === 0) {
                try {
                    var name = reqId.substring(10)
                    var thumb = JSON.parse(data)
                    var updated = Object.assign({}, root._thumbnails)
                    updated[name] = thumb
                    root._thumbnails = updated
                } catch (e) {}
            } else if (reqId === "library_list") {
                try {
                    var libItems = JSON.parse(data)
                    root._libraryAllItems = libItems
                    // Collect unique groups preserving order
                    var seen = {}
                    var groups = ["All"]
                    for (var k = 0; k < libItems.length; ++k) {
                        var g = libItems[k].group || "Other"
                        if (!seen[g]) { seen[g] = true; groups.push(g) }
                    }
                    root._libraryGroups = groups
                    libGroupCombo.model = groups
                    libGroupCombo.currentIndex = 0
                    root._applyLibraryFilter()
                    root._thumbnails = ({})
                    for (var l = 0; l < libItems.length; ++l) {
                        sketch.requestTemplateThumbnail(libItems[l].label, "thumb_lib_" + libItems[l].label)
                    }
                } catch (e) {}
            }
        }
    }

    function _applyLibraryFilter() {
        var group = libGroupCombo.currentText
        var search = libSearchField.text.trim().toLowerCase()
        var all = root._libraryAllItems
        var filtered = []
        for (var i = 0; i < all.length; ++i) {
            if (group !== "All" && all[i].group !== group) continue
            if (search !== "" && all[i].label.toLowerCase().indexOf(search) === -1) continue
            filtered.push(all[i].label)
        }
        libItemList.model = filtered
    }

    AnchoredPicker {
        id: ssPopup
        width: 240
        height: 320
        title: "Salts & Solvents"
        loading: saltsList.count === 0

        GridView {
            id: saltsList
            anchors.fill: parent
            clip: true
            cellWidth: 110
            cellHeight: 48
            model: []
            ScrollBar.vertical: ScrollBar { policy: ScrollBar.AsNeeded }
            delegate: ThumbnailGridItem {
                required property string modelData
                required property int index
                label: modelData
                thumbData: root._thumbnails[modelData] || null
                onClicked: {
                    root.toolSelected("SS_" + modelData)
                    ssPopup.close()
                }
            }
        }
    }

    AnchoredPicker {
        id: fgPopup
        width: 240
        height: 320
        title: "Functional Groups"
        loading: fgList.count === 0

        SegmentedControl {
            id: fgToggleRow
            anchors.top: parent.top
            anchors.left: parent.left
            anchors.topMargin: 4
            anchors.leftMargin: 4
            model: ["Full structure", "Shortcut"]
            currentValue: AppController.fgFullStructure ? "Full structure" : "Shortcut"
            onValueSelected: (value) => { AppController.fgFullStructure = (value === "Full structure") }
        }

        GridView {
            id: fgList
            anchors.top: fgToggleRow.bottom
            anchors.topMargin: 4
            anchors.bottom: parent.bottom
            anchors.left: parent.left
            anchors.right: parent.right
            clip: true
            cellWidth: 110
            cellHeight: 48
            model: []
            ScrollBar.vertical: ScrollBar { policy: ScrollBar.AsNeeded }
            delegate: ThumbnailGridItem {
                required property string modelData
                required property int index
                label: modelData
                thumbData: root._thumbnails[modelData] || null
                onClicked: {
                    root.toolSelected("FG_" + modelData)
                    fgPopup.close()
                }
            }
        }
    }

    // Template Library popup — group filter + item list
    AnchoredPicker {
        id: libPopup
        width: 280
        height: 400
        title: "Template Library"
        loading: libItemList.count === 0

        ColumnLayout {
            anchors.fill: parent
            spacing: 4

            TextField {
                id: libSearchField
                Layout.fillWidth: true
                placeholderText: "Search templates…"
                font.pixelSize: Theme.fontSizeLabel
                onTextChanged: root._applyLibraryFilter()
            }

            ComboBox {
                id: libGroupCombo
                Layout.fillWidth: true
                model: ["All"]
                font.pixelSize: Theme.fontSizeLabel
                onCurrentTextChanged: root._applyLibraryFilter()
            }

            GridView {
                id: libItemList
                Layout.fillWidth: true
                Layout.fillHeight: true
                clip: true
                cellWidth: 110
                cellHeight: 48
                model: []
                ScrollBar.vertical: ScrollBar { policy: ScrollBar.AsNeeded }
                delegate: ThumbnailGridItem {
                    required property string modelData
                    required property int index
                    label: modelData
                    thumbData: root._thumbnails[modelData] || null
                    onClicked: {
                        root.toolSelected("LIB_" + modelData)
                        libPopup.close()
                    }
                }
            }
        }
    }

}
