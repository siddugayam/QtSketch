import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import QtCore
import "../components"
import ".."

Item {
    id: root
    anchors.fill: parent
    required property var win

    property alias settingsDialog: settingsDialog

    // Settings dialog for drawing style and page settings
    TaskDialog {
        id: settingsDialog
        title: "Drawing Settings"
        width: Theme.dialogWidthLarge
        height: Theme.dialogHeightLarge
        standardButtons: Dialog.Apply | Dialog.Cancel

        // SettingsService is a QML_SINGLETON (registered as "SettingsService",
        // no "_cpp" suffix) -- accessed directly as a singleton, never
        // instantiated. See AppController's identical usage elsewhere for
        // the established convention this project already uses.
        property int currentDocId: -1

        // Initialize when opened
        onOpened: {
            currentDocId = DocumentManager.activeDocId
            const doc = DocumentManager.documentFor(currentDocId)
            if (doc) {
                const style = doc.getDrawingStyle()
                const pageSettings = doc.getPageSettings()
                SettingsService.initializeDraftFromJson(style, pageSettings)
            }
        }

        // Apply settings - call target doc via V8Process
        // Uses onApplied (not onAccepted): standardButtons has Dialog.Apply | Dialog.Cancel,
        // and Qt Quick Dialog emits applied() for Apply-role buttons and keeps dialog open.
        onApplied: {
            const doc = DocumentManager.documentFor(currentDocId)
            if (doc) {
                const draftStyle = SettingsService.getDraftDrawingStyleJson()
                const draftPageSettings = SettingsService.getDraftPageSettingsJson()
                doc.applyDrawingStyleWithRescale(draftStyle)
                doc.setPageSettings(draftPageSettings)
                // Initialize draft with the new (now applied) values
                const style = doc.getDrawingStyle()
                const pageSettings = doc.getPageSettings()
                SettingsService.initializeDraftFromJson(style, pageSettings)
            }
        }

        // Cancel - restore original
        onRejected: {
            SettingsService.cancel()
        }

        ColumnLayout {
            anchors.fill: parent
            anchors.margins: Theme.dialogPadding
            spacing: Theme.dialogSpacing

            // Tabs for different setting categories
            TabBar {
                id: settingsTabs
                Layout.fillWidth: true

                TabButton { text: "Drawing Style" }
                TabButton { text: "Page Settings" }
                TabButton { text: "Presets" }
            }

            // Stack layout for tab content
            StackLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                currentIndex: settingsTabs.currentIndex

                // Tab 0: Drawing Style
                ColumnLayout {
                    spacing: Theme.dialogSpacing

                    GroupBox {
                        Layout.fillWidth: true
                        title: "Bond Settings"

                        ColumnLayout {
                            spacing: Theme.dialogSpacing

                            RowLayout {
                                Label { text: "Bond Length:" }
                                TextField {
                                    Layout.fillWidth: true
                                    text: SettingsService.draftBondLength.toFixed(2)
                                    validator: DoubleValidator { bottom: 0.1; top: 10; decimals: 2 }
                                    onAccepted: {
                                        if (acceptableInput) {
                                            SettingsService.draftBondLength = parseFloat(text)
                                        }
                                    }
                                    Keys.onReturnPressed: forceActiveFocus()
                                }
                                Label { text: "model units" }
                            }

                            RowLayout {
                                Label { text: "Line Width:" }
                                TextField {
                                    Layout.fillWidth: true
                                    text: SettingsService.draftLineWidth.toFixed(2)
                                    validator: DoubleValidator { bottom: 0.1; top: 5; decimals: 2 }
                                    onAccepted: {
                                        if (acceptableInput) {
                                            SettingsService.draftLineWidth = parseFloat(text)
                                        }
                                    }
                                }
                                Label { text: "pt" }
                            }
                        }
                    }

                    GroupBox {
                        Layout.fillWidth: true
                        title: "Text Settings"

                        ColumnLayout {
                            spacing: Theme.dialogSpacing

                            RowLayout {
                                Label { text: "Font Family:" }
                                ComboBox {
                                    Layout.fillWidth: true
                                    model: SettingsService.getAvailableFonts()
                                    currentText: SettingsService.draftFontFamily
                                    onCurrentTextChanged: {
                                        SettingsService.draftFontFamily = currentText
                                    }
                                }
                            }

                            RowLayout {
                                Label { text: "Font Size:" }
                                SpinBox {
                                    Layout.fillWidth: true
                                    from: 6
                                    to: 72
                                    value: SettingsService.draftFontSize
                                    onValueChanged: SettingsService.draftFontSize = value
                                }
                                Label { text: "pt" }
                            }
                        }
                    }

                    GroupBox {
                        Layout.fillWidth: true
                        title: "Color Settings"

                        ColumnLayout {
                            spacing: Theme.dialogSpacing

                            RowLayout {
                                CheckBox {
                                    text: "Monochrome Output"
                                    checked: SettingsService.draftMonochrome
                                    onClicked: SettingsService.draftMonochrome = checked
                                }
                            }
                        }
                    }

                    // Preview section
                    GroupBox {
                        Layout.fillWidth: true
                        title: "Live Preview"
                        Layout.fillHeight: true

                        // Preview area
                        PreviewCanvas {
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            primitives: win.activeSketch ? win.activeSketch.primitives : null
                            draftBondLength: SettingsService.draftBondLength
                            draftLineWidth: SettingsService.draftLineWidth
                            draftFontFamily: SettingsService.draftFontFamily
                            draftFontSize: SettingsService.draftFontSize
                            draftMonochrome: SettingsService.draftMonochrome
                            scale: 2.0
                        }
                    }
                }

                // Tab 1: Page Settings
                ColumnLayout {
                    spacing: Theme.dialogSpacing

                    GroupBox {
                        Layout.fillWidth: true
                        title: "Page Size"

                        ColumnLayout {
                            spacing: Theme.dialogSpacing

                            RowLayout {
                                Label { text: "Width:" }
                                TextField {
                                    Layout.fillWidth: true
                                    text: SettingsService.draftPageWidth.toFixed(1)
                                    validator: DoubleValidator { bottom: 10; top: 1000; decimals: 1 }
                                    onAccepted: {
                                        if (acceptableInput) {
                                            SettingsService.draftPageWidth = parseFloat(text)
                                        }
                                    }
                                }
                                Label { text: "mm" }
                            }

                            RowLayout {
                                Label { text: "Height:" }
                                TextField {
                                    Layout.fillWidth: true
                                    text: SettingsService.draftPageHeight.toFixed(1)
                                    validator: DoubleValidator { bottom: 10; top: 1000; decimals: 1 }
                                    onAccepted: {
                                        if (acceptableInput) {
                                            SettingsService.draftPageHeight = parseFloat(text)
                                        }
                                    }
                                }
                                Label { text: "mm" }
                            }
                        }
                    }

                    GroupBox {
                        Layout.fillWidth: true
                        title: "Margins"

                        ColumnLayout {
                            spacing: Theme.dialogSpacing

                            RowLayout {
                                Label { text: "Left:" }
                                TextField {
                                    width: 100
                                    text: SettingsService.draftMarginLeft.toFixed(1)
                                    validator: DoubleValidator { bottom: 0; top: 100; decimals: 1 }
                                    onAccepted: {
                                        if (acceptableInput) {
                                            SettingsService.draftMarginLeft = parseFloat(text)
                                        }
                                    }
                                }
                                Label { text: "mm" }

                                Label { text: "Right:" }
                                TextField {
                                    width: 100
                                    text: SettingsService.draftMarginRight.toFixed(1)
                                    validator: DoubleValidator { bottom: 0; top: 100; decimals: 1 }
                                    onAccepted: {
                                        if (acceptableInput) {
                                            SettingsService.draftMarginRight = parseFloat(text)
                                        }
                                    }
                                }
                                Label { text: "mm" }
                            }

                            RowLayout {
                                Label { text: "Top:" }
                                TextField {
                                    width: 100
                                    text: SettingsService.draftMarginTop.toFixed(1)
                                    validator: DoubleValidator { bottom: 0; top: 100; decimals: 1 }
                                    onAccepted: {
                                        if (acceptableInput) {
                                            SettingsService.draftMarginTop = parseFloat(text)
                                        }
                                    }
                                }
                                Label { text: "mm" }

                                Label { text: "Bottom:" }
                                TextField {
                                    width: 100
                                    text: SettingsService.draftMarginBottom.toFixed(1)
                                    validator: DoubleValidator { bottom: 0; top: 100; decimals: 1 }
                                    onAccepted: {
                                        if (acceptableInput) {
                                            SettingsService.draftMarginBottom = parseFloat(text)
                                        }
                                    }
                                }
                                Label { text: "mm" }
                            }
                        }
                    }
                }

                // Tab 2: Presets
                ColumnLayout {
                    spacing: Theme.dialogSpacing

                    GroupBox {
                        Layout.fillWidth: true
                        title: "Presets"

                        ColumnLayout {
                            spacing: Theme.dialogSpacing

                            // Preset list
                            ListView {
                                Layout.fillWidth: true
                                height: 200
                                clip: true
                                model: SettingsService.getPresets()
                                delegate: ItemDelegate {
                                    width: ListView.view.width
                                    contentItem: RowLayout {
                                        Label {
                                            text: modelData.name + (modelData.isBuiltIn ? " (Built-in)" : "")
                                            Layout.fillWidth: true
                                        }
                                        Button {
                                            text: "Delete"
                                            visible: !modelData.isBuiltIn
                                            onClicked: SettingsService.deletePreset(index)
                                        }
                                    }
                                    onClicked: SettingsService.loadPreset(index)
                                }
                            }

                            RowLayout {
                                TextField {
                                    id: presetNameField
                                    Layout.fillWidth: true
                                    placeholderText: "New preset name..."
                                    Keys.onReturnPressed: savePresetButton.clicked()
                                }
                                Button {
                                    id: savePresetButton
                                    text: "Save"
                                    onClicked: {
                                        const name = presetNameField.text.trim()
                                        if (name !== "") {
                                            if (SettingsService.savePreset(name)) {
                                                presetNameField.text = ""
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    GroupBox {
                        Layout.fillWidth: true
                        title: "Actions"

                        RowLayout {
                            Button {
                                text: "Reset to Defaults"
                                onClicked: SettingsService.resetDraftToDefaults()
                            }
                        }
                    }
                }
            }
        }
    }

    MessageDialog {
        id: presetOverwriteDialog
        title: "Overwrite Preset"
        text: "Preset '" + pendingName + "' already exists. Overwrite?"
        buttons: MessageDialog.Yes | MessageDialog.No
        property string pendingName: ""
        onButtonClicked: function(button) {
            if (button === MessageDialog.Yes) {
                SettingsService.savePreset(pendingName, true)
                pendingName = ""
                presetNameField.text = ""
            }
        }
    }

    MessageDialog {
        id: presetErrorDialog
        title: "Preset Error"
        property string errorText: ""
        text: errorText
        buttons: MessageDialog.Ok
    }

    Connections {
        target: SettingsService
        function onPresetOverwritePrompt(name) {
            presetOverwriteDialog.pendingName = name
            presetOverwriteDialog.open()
        }
        function onPresetError(message) {
            presetErrorDialog.errorText = message
            presetErrorDialog.open()
        }
    }
}
