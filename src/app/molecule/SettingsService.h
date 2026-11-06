// src/app/molecule/SettingsService.h
#ifndef SETTINGSSERVICE_H
#define SETTINGSSERVICE_H

#include <QObject>
#include <QList>
#include <QString>
#include <QVariantMap>
#include <QtQml/qqml.h>
#include "DrawingStyle.h"
#include "PageSettings.h"

// Preset structure for storing named style+page settings
struct StylePreset {
    QString name;
    DrawingStyle drawingStyle;
    PageSettings pageSettings;
};

// QML-accessible settings service singleton
// Manages draft settings and presets for the Drawing Settings dialog
class SettingsService : public QObject {
    Q_OBJECT
    QML_ELEMENT
    QML_SINGLETON

public:
    explicit SettingsService(QObject* parent = nullptr);
    ~SettingsService();

    // Initialize draft from document's current style
    Q_INVOKABLE void initializeDraft(const DrawingStyle& style, const PageSettings& pageSettings);

    // Draft properties (accessed from QML)
    Q_PROPERTY(double draftBondLength READ draftBondLength WRITE setDraftBondLength NOTIFY draftChanged)
    Q_PROPERTY(double draftLineWidth READ draftLineWidth WRITE setDraftLineWidth NOTIFY draftChanged)
    Q_PROPERTY(QString draftFontFamily READ draftFontFamily WRITE setDraftFontFamily NOTIFY draftChanged)
    Q_PROPERTY(int draftFontSize READ draftFontSize WRITE setDraftFontSize NOTIFY draftChanged)
    Q_PROPERTY(bool draftMonochrome READ draftMonochrome WRITE setDraftMonochrome NOTIFY draftChanged)
    Q_PROPERTY(double draftPageWidth READ draftPageWidth WRITE setDraftPageWidth NOTIFY draftChanged)
    Q_PROPERTY(double draftPageHeight READ draftPageHeight WRITE setDraftPageHeight NOTIFY draftChanged)
    Q_PROPERTY(double draftMarginLeft READ draftMarginLeft WRITE setDraftMarginLeft NOTIFY draftChanged)
    Q_PROPERTY(double draftMarginRight READ draftMarginRight WRITE setDraftMarginRight NOTIFY draftChanged)
    Q_PROPERTY(double draftMarginTop READ draftMarginTop WRITE setDraftMarginTop NOTIFY draftChanged)
    Q_PROPERTY(double draftMarginBottom READ draftMarginBottom WRITE setDraftMarginBottom NOTIFY draftChanged)

    // Getters
    double draftBondLength() const { return m_draftDrawingStyle.bondLength; }
    double draftLineWidth() const { return m_draftDrawingStyle.lineWidth; }
    QString draftFontFamily() const { return m_draftDrawingStyle.fontFamily; }
    int draftFontSize() const { return m_draftDrawingStyle.fontSize; }
    bool draftMonochrome() const { return m_draftDrawingStyle.monochrome; }
    double draftPageWidth() const { return m_draftPageSettings.widthMm; }
    double draftPageHeight() const { return m_draftPageSettings.heightMm; }
    double draftMarginLeft() const { return m_draftPageSettings.marginLeftMm; }
    double draftMarginRight() const { return m_draftPageSettings.marginRightMm; }
    double draftMarginTop() const { return m_draftPageSettings.marginTopMm; }
    double draftMarginBottom() const { return m_draftPageSettings.marginBottomMm; }

    // Setters
    void setDraftBondLength(double value);
    void setDraftLineWidth(double value);
    void setDraftFontFamily(const QString& value);
    void setDraftFontSize(int value);
    void setDraftMonochrome(bool value);
    void setDraftPageWidth(double value);
    void setDraftPageHeight(double value);
    void setDraftMarginLeft(double value);
    void setDraftMarginRight(double value);
    void setDraftMarginTop(double value);
    void setDraftMarginBottom(double value);

    // Get the current draft as style/page settings objects
    Q_INVOKABLE DrawingStyle getDraftDrawingStyle() const;
    Q_INVOKABLE PageSettings getDraftPageSettings() const;
    
    // Get draft as QJsonObject for easy QML→C++ passing
    Q_INVOKABLE QJsonObject getDraftDrawingStyleJson() const;
    Q_INVOKABLE QJsonObject getDraftPageSettingsJson() const;

    // Preset management
    Q_INVOKABLE QVariantList getPresets() const;
    Q_INVOKABLE bool loadPreset(int index);
    Q_INVOKABLE bool savePreset(const QString& name);
    Q_INVOKABLE bool deletePreset(int index);
    Q_INVOKABLE QStringList getAvailableFonts() const;

    // Apply draft settings to document (returns the new style and page settings)
    // This triggers the rescale if bond length changed
    Q_INVOKABLE void applyToDocument(int docId);

    // Cancel and reset draft
    Q_INVOKABLE void cancel();

    // Check if draft differs from original
    Q_INVOKABLE bool hasChanges() const;

    // Reset draft to defaults
    Q_INVOKABLE void resetDraftToDefaults();

Q_SIGNALS:
    void draftChanged();
    void presetsChanged();
    void applyRequested(int docId);
    void cancelRequested();

private:
    // Load presets from storage
    void loadPresetsFromStorage();

    // Save presets to storage
    void savePresetsToStorage();

    // Get presets storage path
    QString getPresetsPath() const;

    // Draft settings (not applied until Apply is pressed)
    DrawingStyle m_draftDrawingStyle;
    PageSettings m_draftPageSettings;

    // Original settings (for Cancel)
    DrawingStyle m_originalDrawingStyle;
    PageSettings m_originalPageSettings;

    // Loaded presets
    QList<StylePreset> m_presets;

    // Track if presets are loaded
    bool m_presetsLoaded = false;
};

#endif // SETTINGSSERVICE_H
