// tests/settings_service_test.cpp
// Tests for SettingsService class
#include <cstdio>
#include <QGuiApplication>
#include <QSignalSpy>
#include <QStandardPaths>
#include <QList>
#include "app/molecule/SettingsService.h"
#include "app/molecule/DrawingStyle.h"
#include "app/molecule/PageSettings.h"

static int g_pass = 0, g_fail = 0;
#define CHECK(cond, name) do { \
    if (cond) { ++g_pass; std::printf("[PASS] %s\n", name); } \
    else { ++g_fail; std::printf("[FAIL] %s (line %d)\n", name, __LINE__); } \
} while (0)

int main() {
    // SettingsService links Qt6::Gui transitively; use QGuiApplication (not
    // QCoreApplication) -- confirmed real fail-fast crash (0xC0000602) with
    // QCoreApplication before this fix.
    int argc = 0;
    char** argv = nullptr;
    QGuiApplication app(argc, argv);
    
    // Enable test mode for QStandardPaths to avoid writing to real user data
    QStandardPaths::setTestModeEnabled(true);

    std::printf("--- SettingsService Tests ---\n\n");

    // Create SettingsService instance
    SettingsService service;

    // Test 1: initializeDraft correctly seeds draft and original
    std::printf("--- Test 1: initializeDraft ---\n");
    DrawingStyle style1(1.5, 0.6, "Arial", 10, true);
    PageSettings page1(210.0, 297.0, 5.0, 5.0, 5.0, 5.0);
    service.initializeDraft(style1, page1);
    
    CHECK(service.draftBondLength() == 1.5, "draft bond length initialized");
    CHECK(service.draftLineWidth() == 0.6, "draft line width initialized");
    CHECK(service.draftFontFamily() == "Arial", "draft font family initialized");
    CHECK(service.draftFontSize() == 10, "draft font size initialized");
    CHECK(service.draftMonochrome() == true, "draft monochrome initialized");
    CHECK(service.draftPageWidth() == 210.0, "draft page width initialized");
    CHECK(service.draftPageHeight() == 297.0, "draft page height initialized");
    CHECK(service.draftMarginLeft() == 5.0, "draft margin left initialized");
    CHECK(service.draftMarginRight() == 5.0, "draft margin right initialized");
    CHECK(service.draftMarginTop() == 5.0, "draft margin top initialized");
    CHECK(service.draftMarginBottom() == 5.0, "draft margin bottom initialized");

    // Test 2: hasChanges returns false immediately after initializeDraft
    std::printf("--- Test 2: hasChanges after initialize ---\n");
    CHECK(!service.hasChanges(), "hasChanges returns false after initializeDraft");

    // Test 3: Setters update fields and emit draftChanged
    std::printf("--- Test 3: Setters emit draftChanged ---\n");
    QSignalSpy draftChangedSpy(&service, &SettingsService::draftChanged);
    
    service.setDraftBondLength(2.0);
    CHECK(draftChangedSpy.count() == 1, "setDraftBondLength emits draftChanged");
    CHECK(service.draftBondLength() == 2.0, "draft bond length updated");
    CHECK(service.hasChanges(), "hasChanges returns true after setter call");
    
    service.setDraftLineWidth(0.8);
    CHECK(draftChangedSpy.count() == 2, "setDraftLineWidth emits draftChanged");
    CHECK(service.draftLineWidth() == 0.8, "draft line width updated");
    
    service.setDraftFontFamily("Times New Roman");
    CHECK(draftChangedSpy.count() == 3, "setDraftFontFamily emits draftChanged");
    CHECK(service.draftFontFamily() == "Times New Roman", "draft font family updated");
    
    service.setDraftFontSize(12);
    CHECK(draftChangedSpy.count() == 4, "setDraftFontSize emits draftChanged");
    CHECK(service.draftFontSize() == 12, "draft font size updated");
    
    service.setDraftMonochrome(false);
    CHECK(draftChangedSpy.count() == 5, "setDraftMonochrome emits draftChanged");
    CHECK(service.draftMonochrome() == false, "draft monochrome updated");
    
    service.setDraftPageWidth(300.0);
    CHECK(draftChangedSpy.count() == 6, "setDraftPageWidth emits draftChanged");
    CHECK(service.draftPageWidth() == 300.0, "draft page width updated");
    
    service.setDraftPageHeight(400.0);
    CHECK(draftChangedSpy.count() == 7, "setDraftPageHeight emits draftChanged");
    CHECK(service.draftPageHeight() == 400.0, "draft page height updated");
    
    service.setDraftMarginLeft(10.0);
    CHECK(draftChangedSpy.count() == 8, "setDraftMarginLeft emits draftChanged");
    CHECK(service.draftMarginLeft() == 10.0, "draft margin left updated");
    
    service.setDraftMarginRight(15.0);
    CHECK(draftChangedSpy.count() == 9, "setDraftMarginRight emits draftChanged");
    
    service.setDraftMarginTop(20.0);
    CHECK(draftChangedSpy.count() == 10, "setDraftMarginTop emits draftChanged");
    
    service.setDraftMarginBottom(25.0);
    CHECK(draftChangedSpy.count() == 11, "setDraftMarginBottom emits draftChanged");

    // Test 4: cancel resets draft to original and emits cancelRequested
    std::printf("--- Test 4: cancel resets to original ---\n");
    QSignalSpy cancelRequestedSpy(&service, &SettingsService::cancelRequested);
    service.cancel();
    CHECK(draftChangedSpy.count() == 12, "cancel emits draftChanged");
    CHECK(cancelRequestedSpy.count() == 1, "cancel emits cancelRequested");
    CHECK(service.draftBondLength() == 1.5, "draft bond length reset to original");
    CHECK(service.draftLineWidth() == 0.6, "draft line width reset to original");
    CHECK(service.draftFontFamily() == "Arial", "draft font family reset to original");
    CHECK(service.draftFontSize() == 10, "draft font size reset to original");
    CHECK(service.draftMonochrome() == true, "draft monochrome reset to original");
    CHECK(service.draftPageWidth() == 210.0, "draft page width reset to original");
    CHECK(service.draftPageHeight() == 297.0, "draft page height reset to original");
    CHECK(service.draftMarginLeft() == 5.0, "draft margin left reset to original");
    CHECK(!service.hasChanges(), "hasChanges returns false after cancel");

    // Test 5: resetDraftToDefaults
    std::printf("--- Test 5: resetDraftToDefaults ---\n");
    service.resetDraftToDefaults();
    DrawingStyle defaults = DrawingStyle::defaults();
    PageSettings pageDefaults = PageSettings::defaults();
    CHECK(service.draftBondLength() == defaults.bondLength, "draft bond length reset to default");
    CHECK(service.draftLineWidth() == defaults.lineWidth, "draft line width reset to default");
    CHECK(service.draftFontFamily() == defaults.fontFamily, "draft font family reset to default");
    CHECK(service.draftFontSize() == defaults.fontSize, "draft font size reset to default");
    CHECK(service.draftMonochrome() == defaults.monochrome, "draft monochrome reset to default");
    CHECK(service.draftPageWidth() == pageDefaults.widthMm, "draft page width reset to default");
    CHECK(service.draftPageHeight() == pageDefaults.heightMm, "draft page height reset to default");
    CHECK(service.draftMarginLeft() == pageDefaults.marginLeftMm, "draft margin left reset to default");
    CHECK(service.draftMarginRight() == pageDefaults.marginRightMm, "draft margin right reset to default");
    CHECK(service.draftMarginTop() == pageDefaults.marginTopMm, "draft margin top reset to default");
    CHECK(service.draftMarginBottom() == pageDefaults.marginBottomMm, "draft margin bottom reset to default");

    // Test 6: applyToDocument emits applyRequested with docId
    std::printf("--- Test 6: applyToDocument emits applyRequested ---\n");
    QSignalSpy applyRequestedSpy(&service, &SettingsService::applyRequested);
    service.applyToDocument(42);
    CHECK(applyRequestedSpy.count() == 1, "applyToDocument emits applyRequested");
    QList<QVariant> applyArgs = applyRequestedSpy.takeFirst();
    CHECK(applyArgs.size() == 1, "applyRequested has one argument");
    CHECK(applyArgs[0].toInt() == 42, "applyRequested argument is the docId");

    // Test 7: Preset round-trip
    std::printf("--- Test 7: Preset round-trip ---\n");
    // Initialize with known values
    service.initializeDraft(DrawingStyle::defaults(), PageSettings::defaults());
    
    // Get initial presets (should include built-in ones)
    QVariantList presetsBefore = service.getPresets();
    int initialCount = presetsBefore.size();
    CHECK(initialCount > 0, "initial presets list is not empty");
    
    // Save a new preset
    service.setDraftBondLength(3.0);
    bool saved = service.savePreset("Test Preset");
    CHECK(saved, "savePreset returns true");
    
    // Check preset is in the list
    QVariantList presetsAfterSave = service.getPresets();
    CHECK(presetsAfterSave.size() == initialCount + 1, "preset added to list");
    
    // Load the preset
    bool loaded = service.loadPreset(presetsAfterSave.size() - 1);  // Last one should be our test preset
    CHECK(loaded, "loadPreset returns true");
    CHECK(service.draftBondLength() == 3.0, "loaded preset has correct bond length");
    
    // Delete the preset
    bool deleted = service.deletePreset(presetsAfterSave.size() - 1);
    CHECK(deleted, "deletePreset returns true");
    
    QVariantList presetsAfterDelete = service.getPresets();
    CHECK(presetsAfterDelete.size() == initialCount, "preset removed from list");
    // Verify the preset is gone by trying to load it (should fail or have different value)
    service.initializeDraft(DrawingStyle::defaults(), PageSettings::defaults());
    bool shouldFail = !service.loadPreset(presetsAfterSave.size() - 1);  // Index should be out of bounds now
    CHECK(shouldFail, "loadPreset fails for deleted preset index");

    // Test 8: getAvailableFonts returns non-empty list
    std::printf("--- Test 8: getAvailableFonts ---\n");
    QStringList fonts = service.getAvailableFonts();
    CHECK(!fonts.isEmpty(), "getAvailableFonts returns non-empty list");

    // Test 9: initializeDraftFromJson converts QJsonObject to draft properties
    std::printf("--- Test 9: initializeDraftFromJson ---\n");
    QJsonObject styleJson;
    styleJson["bondLength"] = 2.5;
    styleJson["lineWidth"] = 0.75;
    styleJson["fontFamily"] = "Helvetica";
    styleJson["fontSize"] = 14;
    styleJson["monochrome"] = false;
    styleJson["modelToMmScale"] = 3.5;

    QJsonObject pageJson;
    pageJson["widthMm"] = 150.0;
    pageJson["heightMm"] = 250.0;
    pageJson["marginLeftMm"] = 12.0;
    pageJson["marginRightMm"] = 13.0;
    pageJson["marginTopMm"] = 14.0;
    pageJson["marginBottomMm"] = 15.0;

    service.initializeDraftFromJson(styleJson, pageJson);
    CHECK(service.draftBondLength() == 2.5, "draft bond length from JSON");
    CHECK(service.draftLineWidth() == 0.75, "draft line width from JSON");
    CHECK(service.draftFontFamily() == "Helvetica", "draft font family from JSON");
    CHECK(service.draftFontSize() == 14, "draft font size from JSON");
    CHECK(service.draftMonochrome() == false, "draft monochrome from JSON");
    CHECK(service.draftPageWidth() == 150.0, "draft page width from JSON");
    CHECK(service.draftPageHeight() == 250.0, "draft page height from JSON");
    CHECK(service.draftMarginLeft() == 12.0, "draft margin left from JSON");
    CHECK(service.draftMarginRight() == 13.0, "draft margin right from JSON");
    CHECK(service.draftMarginTop() == 14.0, "draft margin top from JSON");
    CHECK(service.draftMarginBottom() == 15.0, "draft margin bottom from JSON");

    // Test 10: Preset collision and built-in protection
    std::printf("--- Test 10: Preset protection ---\n");
    QSignalSpy errorSpy(&service, &SettingsService::presetError);
    QSignalSpy promptSpy(&service, &SettingsService::presetOverwritePrompt);

    // Built-in deletion rejected
    CHECK(!service.deletePreset(0), "deletePreset rejects built-in index 0");
    CHECK(errorSpy.count() == 1, "deletePreset emits presetError for built-in");

    // Built-in overwrite rejected
    CHECK(!service.savePreset("Default"), "savePreset rejects overwriting built-in Default");
    CHECK(errorSpy.count() == 2, "savePreset emits presetError for built-in overwrite attempt");

    // User preset collision prompt
    CHECK(service.savePreset("CustomPreset"), "save new CustomPreset succeeds");
    CHECK(!service.savePreset("CustomPreset", false), "save colliding CustomPreset without allowOverwrite rejected");
    CHECK(promptSpy.count() == 1, "colliding savePreset emits presetOverwritePrompt");
    CHECK(promptSpy.first().at(0).toString() == "CustomPreset", "prompt argument is preset name");

    // User preset overwrite confirmed
    service.setDraftBondLength(4.0);
    CHECK(service.savePreset("CustomPreset", true), "save colliding CustomPreset with allowOverwrite=true succeeds");
    CHECK(service.loadPreset(service.getPresets().size() - 1), "load updated CustomPreset succeeds");
    CHECK(service.draftBondLength() == 4.0, "updated preset has new bond length");

    // Clean up created user preset
    CHECK(service.deletePreset(service.getPresets().size() - 1), "delete custom preset succeeds");

    std::printf("\nSummary: %d passed, %d failed.\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
