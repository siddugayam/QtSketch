import QtQuick
import QtQuick.Layouts
import QtQuick.Templates as T
import Sketch.App

T.AbstractButton {
    id: control
    property string shortcutHint: ""
    property string iconSource: ""

    signal triggered()

    implicitWidth: Math.max(220, contentRow.implicitWidth + 24)
    implicitHeight: 28

    Accessible.role: Accessible.MenuItem
    Keys.onReturnPressed: clicked()
    Keys.onSpacePressed: clicked()

    onClicked: {
        triggered()
        // Find the popup to close it
        var p = parent
        while (p && typeof p.close === "undefined") {
            p = p.parent
        }
        if (p) p.close()
    }

    background: Rectangle {
        radius: 3
        color: control.down ? Theme.selected : (control.hovered ? Theme.hover : "transparent")
    }

    contentItem: Item {
        anchors.fill: parent
        anchors.leftMargin: 8
        anchors.rightMargin: 8

        RowLayout {
            id: contentRow
            anchors.fill: parent
            spacing: 8

            Image {
                visible: control.iconSource !== ""
                source: control.iconSource !== "" ? "../icons/" + control.iconSource : ""
                Layout.preferredWidth: 16
                Layout.preferredHeight: 16
                sourceSize: Qt.size(16, 16)
            }

            Text {
                text: control.text
                color: Theme.textPrimary
                font { pixelSize: Theme.fontSizeBody; family: Theme.fontFamily }
                Layout.fillWidth: true
            }

            Text {
                visible: control.shortcutHint !== ""
                text: Theme.formatShortcut(control.shortcutHint)
                color: Theme.textSecondary
                font { pixelSize: Theme.fontSizeBody; family: Theme.fontFamily }
                horizontalAlignment: Text.AlignRight
            }
        }
    }

    FocusRing {
        anchors.fill: parent
        visible: control.activeFocus
        radius: 3
    }
}
