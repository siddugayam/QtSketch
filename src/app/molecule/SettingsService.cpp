// src/app/molecule/SettingsService.cpp
#include "SettingsService.h"
#include "DocumentManager.h"
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QStandardPaths>
#include <QFontDatabase>
#include <QDebug>
#include <QDir>

SettingsService::SettingsService(QObject* parent) : QObject(parent) {
    // Load presets from storage
    loadPresetsFromStorage();
    
    // Initialize with defaults
    m_draftDrawingStyle = DrawingStyle::defaults();
    m_draftPageSettings = PageSettings::defaults();
    m_originalDrawingStyle = m_draftDrawingStyle;
    m_originalPageSettings = m_draftPageSettings;
}

SettingsService::~SettingsService() {
    savePresetsToStorage();
}

void SettingsService::initializeDraft(const DrawingStyle& style, const PageSettings& pageSettings) {
    m_originalDrawingStyle = style;
    m_originalPageSettings = pageSettings;
    m_draftDrawingStyle = style;
    m_draftPageSettings = pageSettings;
    emit draftChanged();
}

// Setters with emit
void SettingsService::setDraftBondLength(double value) {
    if (m_draftDrawingStyle.bondLength != value) {
        m_draftDrawingStyle.bondLength = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftLineWidth(double value) {
    if (m_draftDrawingStyle.lineWidth != value) {
        m_draftDrawingStyle.lineWidth = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftFontFamily(const QString& value) {
    if (m_draftDrawingStyle.fontFamily != value) {
        m_draftDrawingStyle.fontFamily = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftFontSize(int value) {
    if (m_draftDrawingStyle.fontSize != value) {
        m_draftDrawingStyle.fontSize = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftMonochrome(bool value) {
    if (m_draftDrawingStyle.monochrome != value) {
        m_draftDrawingStyle.monochrome = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftPageWidth(double value) {
    if (m_draftPageSettings.widthMm != value) {
        m_draftPageSettings.widthMm = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftPageHeight(double value) {
    if (m_draftPageSettings.heightMm != value) {
        m_draftPageSettings.heightMm = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftMarginLeft(double value) {
    if (m_draftPageSettings.marginLeftMm != value) {
        m_draftPageSettings.marginLeftMm = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftMarginRight(double value) {
    if (m_draftPageSettings.marginRightMm != value) {
        m_draftPageSettings.marginRightMm = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftMarginTop(double value) {
    if (m_draftPageSettings.marginTopMm != value) {
        m_draftPageSettings.marginTopMm = value;
        emit draftChanged();
    }
}

void SettingsService::setDraftMarginBottom(double value) {
    if (m_draftPageSettings.marginBottomMm != value) {
        m_draftPageSettings.marginBottomMm = value;
        emit draftChanged();
    }
}

DrawingStyle SettingsService::getDraftDrawingStyle() const {
    return m_draftDrawingStyle;
}

PageSettings SettingsService::getDraftPageSettings() const {
    return m_draftPageSettings;
}

QJsonObject SettingsService::getDraftDrawingStyleJson() const {
    QJsonObject obj;
    obj["bondLength"] = m_draftDrawingStyle.bondLength;
    obj["lineWidth"] = m_draftDrawingStyle.lineWidth;
    obj["fontFamily"] = m_draftDrawingStyle.fontFamily;
    obj["fontSize"] = m_draftDrawingStyle.fontSize;
    obj["monochrome"] = m_draftDrawingStyle.monochrome;
    obj["modelToMmScale"] = m_draftDrawingStyle.modelToMmScale;
    return obj;
}

QJsonObject SettingsService::getDraftPageSettingsJson() const {
    QJsonObject obj;
    obj["widthMm"] = m_draftPageSettings.widthMm;
    obj["heightMm"] = m_draftPageSettings.heightMm;
    obj["marginLeftMm"] = m_draftPageSettings.marginLeftMm;
    obj["marginRightMm"] = m_draftPageSettings.marginRightMm;
    obj["marginTopMm"] = m_draftPageSettings.marginTopMm;
    obj["marginBottomMm"] = m_draftPageSettings.marginBottomMm;
    return obj;
}

QVariantList SettingsService::getPresets() const {
    QVariantList result;
    for (const StylePreset& preset : m_presets) {
        QVariantMap map;
        map["name"] = preset.name;
        map["bondLength"] = preset.drawingStyle.bondLength;
        map["lineWidth"] = preset.drawingStyle.lineWidth;
        map["fontFamily"] = preset.drawingStyle.fontFamily;
        map["fontSize"] = preset.drawingStyle.fontSize;
        map["monochrome"] = preset.drawingStyle.monochrome;
        map["pageWidth"] = preset.pageSettings.widthMm;
        map["pageHeight"] = preset.pageSettings.heightMm;
        map["marginLeft"] = preset.pageSettings.marginLeftMm;
        map["marginRight"] = preset.pageSettings.marginRightMm;
        map["marginTop"] = preset.pageSettings.marginTopMm;
        map["marginBottom"] = preset.pageSettings.marginBottomMm;
        result.append(map);
    }
    return result;
}

bool SettingsService::loadPreset(int index) {
    if (index < 0 || index >= m_presets.size()) {
        return false;
    }
    
    const StylePreset& preset = m_presets[index];
    m_draftDrawingStyle = preset.drawingStyle;
    m_draftPageSettings = preset.pageSettings;
    emit draftChanged();
    return true;
}

bool SettingsService::savePreset(const QString& name) {
    if (name.isEmpty()) {
        return false;
    }
    
    // Check if preset with this name already exists
    for (int i = 0; i < m_presets.size(); ++i) {
        if (m_presets[i].name == name) {
            // Update existing preset
            m_presets[i].drawingStyle = m_draftDrawingStyle;
            m_presets[i].pageSettings = m_draftPageSettings;
            savePresetsToStorage();
            emit presetsChanged();
            return true;
        }
    }
    
    // Add new preset
    StylePreset preset;
    preset.name = name;
    preset.drawingStyle = m_draftDrawingStyle;
    preset.pageSettings = m_draftPageSettings;
    m_presets.append(preset);
    savePresetsToStorage();
    emit presetsChanged();
    return true;
}

bool SettingsService::deletePreset(int index) {
    if (index < 0 || index >= m_presets.size()) {
        return false;
    }
    
    m_presets.removeAt(index);
    savePresetsToStorage();
    emit presetsChanged();
    return true;
}

QStringList SettingsService::getAvailableFonts() const {
    QStringList fonts;
    QFontDatabase db;
    fonts = db.families();
    return fonts;
}

void SettingsService::applyToDocument(int docId) {
    // Get the document's V8Process and apply the draft settings
    // This is called from QML when Apply is clicked
    emit applyRequested(docId);
    
    // Also emit that apply was completed
    // The actual application happens in QML via the signal
}

void SettingsService::cancel() {
    m_draftDrawingStyle = m_originalDrawingStyle;
    m_draftPageSettings = m_originalPageSettings;
    emit draftChanged();
    emit cancelRequested();
}

bool SettingsService::hasChanges() const {
    return m_draftDrawingStyle != m_originalDrawingStyle ||
           m_draftPageSettings != m_originalPageSettings;
}

void SettingsService::resetDraftToDefaults() {
    m_draftDrawingStyle = DrawingStyle::defaults();
    m_draftPageSettings = PageSettings::defaults();
    emit draftChanged();
}

QString SettingsService::getPresetsPath() const {
    QString appData = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);
    return appData + "/presets/drawing_style_presets.json";
}

void SettingsService::loadPresetsFromStorage() {
    if (m_presetsLoaded) return;
    
    QString path = getPresetsPath();
    QFile file(path);
    
    if (!file.open(QIODevice::ReadOnly)) {
        // No file yet, use defaults
        // Add some built-in presets
        StylePreset defaultPreset;
        defaultPreset.name = "Default";
        defaultPreset.drawingStyle = DrawingStyle::defaults();
        defaultPreset.pageSettings = PageSettings::defaults();
        m_presets.append(defaultPreset);
        
        StylePreset publicationPreset;
        publicationPreset.name = "Publication";
        publicationPreset.drawingStyle = DrawingStyle(1.5, 0.6, "Arial", 10, true);
        publicationPreset.pageSettings = PageSettings(210, 297, 10, 10, 10, 10);
        m_presets.append(publicationPreset);
        
        StylePreset presentationPreset;
        presentationPreset.name = "Presentation";
        presentationPreset.drawingStyle = DrawingStyle(2.0, 1.0, "Arial", 14, false);
        presentationPreset.pageSettings = PageSettings(210, 297, 20, 20, 20, 20);
        m_presets.append(presentationPreset);
        
        m_presetsLoaded = true;
        return;
    }
    
    QByteArray data = file.readAll();
    file.close();
    
    QJsonDocument doc = QJsonDocument::fromJson(data);
    if (!doc.isArray()) {
        qWarning() << "SettingsService: Invalid presets format, using defaults";
        // Use defaults
        StylePreset defaultPreset;
        defaultPreset.name = "Default";
        defaultPreset.drawingStyle = DrawingStyle::defaults();
        defaultPreset.pageSettings = PageSettings::defaults();
        m_presets.append(defaultPreset);
        m_presetsLoaded = true;
        return;
    }
    
    QJsonArray array = doc.array();
    for (const QJsonValue& value : array) {
        if (value.isObject()) {
            QJsonObject obj = value.toObject();
            StylePreset preset;
            preset.name = obj["name"].toString();
            
            QJsonObject styleObj = obj["drawingStyle"].toObject();
            preset.drawingStyle.bondLength = styleObj["bondLength"].toDouble(1.5);
            preset.drawingStyle.lineWidth = styleObj["lineWidth"].toDouble(0.6);
            preset.drawingStyle.fontFamily = styleObj["fontFamily"].toString("Arial");
            preset.drawingStyle.fontSize = styleObj["fontSize"].toInt(10);
            preset.drawingStyle.monochrome = styleObj["monochrome"].toBool(true);
            preset.drawingStyle.modelToMmScale = styleObj["modelToMmScale"].toDouble(5.08 / 1.5);
            
            QJsonObject pageObj = obj["pageSettings"].toObject();
            preset.pageSettings.widthMm = pageObj["widthMm"].toDouble(210);
            preset.pageSettings.heightMm = pageObj["heightMm"].toDouble(297);
            preset.pageSettings.marginLeftMm = pageObj["marginLeftMm"].toDouble(0);
            preset.pageSettings.marginRightMm = pageObj["marginRightMm"].toDouble(0);
            preset.pageSettings.marginTopMm = pageObj["marginTopMm"].toDouble(0);
            preset.pageSettings.marginBottomMm = pageObj["marginBottomMm"].toDouble(0);
            
            m_presets.append(preset);
        }
    }
    
    m_presetsLoaded = true;
}

void SettingsService::savePresetsToStorage() {
    if (!m_presetsLoaded) return;
    
    QString path = getPresetsPath();
    QDir().mkpath(QFileInfo(path).dir().path());
    
    QJsonArray array;
    for (const StylePreset& preset : m_presets) {
        QJsonObject obj;
        obj["name"] = preset.name;
        
        QJsonObject styleObj;
        styleObj["bondLength"] = preset.drawingStyle.bondLength;
        styleObj["lineWidth"] = preset.drawingStyle.lineWidth;
        styleObj["fontFamily"] = preset.drawingStyle.fontFamily;
        styleObj["fontSize"] = preset.drawingStyle.fontSize;
        styleObj["monochrome"] = preset.drawingStyle.monochrome;
        styleObj["modelToMmScale"] = preset.drawingStyle.modelToMmScale;
        obj["drawingStyle"] = styleObj;
        
        QJsonObject pageObj;
        pageObj["widthMm"] = preset.pageSettings.widthMm;
        pageObj["heightMm"] = preset.pageSettings.heightMm;
        pageObj["marginLeftMm"] = preset.pageSettings.marginLeftMm;
        pageObj["marginRightMm"] = preset.pageSettings.marginRightMm;
        pageObj["marginTopMm"] = preset.pageSettings.marginTopMm;
        pageObj["marginBottomMm"] = preset.pageSettings.marginBottomMm;
        obj["pageSettings"] = pageObj;
        
        array.append(obj);
    }
    
    QJsonDocument doc(array);
    QFile file(path);
    if (file.open(QIODevice::WriteOnly)) {
        file.write(doc.toJson());
        file.close();
    } else {
        qWarning() << "SettingsService: Failed to save presets to" << path;
    }
}
