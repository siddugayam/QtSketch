// src/app/molecule/ProjectFile.cpp
#include "ProjectFile.h"
#include "EditableMolecule.h"
#include <QFile>
#include <QSaveFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QDebug>
#include <QRegularExpression>

namespace {
    // Keys for JSON serialization
    const QString kFormatVersion = QStringLiteral("formatVersion");
    const QString kChemistry = QStringLiteral("chemistry");
    const QString kMolfile = QStringLiteral("molfile");
    const QString kDrawingStyle = QStringLiteral("drawingStyle");
    const QString kPageSettings = QStringLiteral("pageSettings");
    const QString kModelToMmScale = QStringLiteral("modelToMmScale");
    const QString kEzIntent = QStringLiteral("ezIntent");
    const QString kSgroups = QStringLiteral("sgroups");
    const QString kAnnotations = QStringLiteral("annotations");
    const QString kImages = QStringLiteral("images");
    const QString kRxnArrows = QStringLiteral("rxnArrows");
    const QString kRxnPluses = QStringLiteral("rxnPluses");
    const QString kMultitailArrows = QStringLiteral("multitailArrows");
    const QString kRGroups = QStringLiteral("rGroups");
    const QString kBrackets = QStringLiteral("brackets");
    const QString kName = QStringLiteral("name");
    const QString kStereoFlags = QStringLiteral("stereoFlags");
    
    // DrawingStyle keys
    const QString kBondLength = QStringLiteral("bondLength");
    const QString kLineWidth = QStringLiteral("lineWidth");
    const QString kFontFamily = QStringLiteral("fontFamily");
    const QString kFontSize = QStringLiteral("fontSize");
    const QString kMonochrome = QStringLiteral("monochrome");
    
    // PageSettings keys
    const QString kWidthMm = QStringLiteral("widthMm");
    const QString kHeightMm = QStringLiteral("heightMm");
    const QString kMarginLeftMm = QStringLiteral("marginLeftMm");
    const QString kMarginRightMm = QStringLiteral("marginRightMm");
    const QString kMarginTopMm = QStringLiteral("marginTopMm");
    const QString kMarginBottomMm = QStringLiteral("marginBottomMm");
    
    // Annotation keys
    const QString kType = QStringLiteral("type");
    const QString kId = QStringLiteral("id");
    const QString kX = QStringLiteral("x");
    const QString kY = QStringLiteral("y");
    const QString kContent = QStringLiteral("content");
    const QString kBold = QStringLiteral("bold");
    const QString kItalic = QStringLiteral("italic");
    const QString kW = QStringLiteral("w");
    const QString kH = QStringLiteral("h");
    const QString kData = QStringLiteral("data");
    const QString kPoints = QStringLiteral("points");
    const QString kMode = QStringLiteral("mode");
    const QString kAbove = QStringLiteral("above");
    const QString kBelow = QStringLiteral("below");
    const QString kNumber = QStringLiteral("number");
    const QString kRange = QStringLiteral("range");
    const QString kResth = QStringLiteral("resth");
    const QString kIfthen = QStringLiteral("ifthen");
    const QString kFragments = QStringLiteral("fragments");
    const QString kAtoms = QStringLiteral("atoms");
    const QString kExpanded = QStringLiteral("expanded");
    const QString kLabel = QStringLiteral("label");
    const QString kMinX = QStringLiteral("minX");
    const QString kMinY = QStringLiteral("minY");
    const QString kMaxX = QStringLiteral("maxX");
    const QString kMaxY = QStringLiteral("maxY");
    
    // Check if a string looks like a valid molfile (has atom/bond counts)
    bool looksLikeMolfile(const QString& data) {
        // Simple heuristic: molfiles typically have "  0  0  0  0  0  0" pattern
        QRegularExpression molfilePattern(R"(^\s*\d+\s+\d+\s+\d+)");
        return molfilePattern.match(data).hasMatch();
    }
}

bool ProjectFile::saveProject(const QString& filePath,
                              const EditableMolecule& molecule,
                              const DrawingStyle& drawingStyle,
                              const PageSettings& pageSettings) {
    QJsonObject json = toJson(molecule, drawingStyle, pageSettings);
    QByteArray bytes = QJsonDocument(json).toJson(QJsonDocument::Indented);

    // QSaveFile writes to its own temp file and only replaces the real
    // destination in commit() once the full write is confirmed -- unlike the
    // old manual temp-file-then-rename dance, it never deletes the existing
    // file up front, so a write/commit failure always leaves the previous
    // project exactly as it was.
    QSaveFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qWarning() << "ProjectFile::saveProject: Failed to open for writing:" << filePath;
        return false;
    }

    const qint64 written = file.write(bytes);
    if (written != bytes.size()) {
        qWarning() << "ProjectFile::saveProject: Partial or failed write to:" << filePath
                   << "(" << written << "of" << bytes.size() << "bytes)";
        file.cancelWriting();
        return false;
    }

    if (!file.commit()) {
        qWarning() << "ProjectFile::saveProject: Failed to commit:" << filePath
                   << "-" << file.errorString();
        return false;
    }

    return true;
}

QString ProjectFile::loadProject(const QString& filePath,
                                  EditableMolecule& molecule,
                                  DrawingStyle& drawingStyle,
                                  PageSettings& pageSettings) {
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return QStringLiteral("Failed to open file: ") + filePath;
    }
    
    QByteArray data = file.readAll();
    file.close();
    
    QJsonParseError parseError;
    QJsonDocument doc = QJsonDocument::fromJson(data, &parseError);
    
    if (parseError.error != QJsonParseError::NoError) {
        return QStringLiteral("Invalid JSON in file: ") + parseError.errorString();
    }
    
    if (!doc.isObject()) {
        return QStringLiteral("File does not contain a JSON object");
    }
    
    QJsonObject json = doc.object();

    // Validate the full schema before touching any live state -- a
    // structurally invalid file (missing version, missing chemistry) must
    // never reach fromJson's mutations at all.
    QString validationError;
    if (!validateProjectJson(json, validationError)) {
        return validationError;
    }

    if (!fromJson(json, molecule, drawingStyle, pageSettings)) {
        return QStringLiteral("Failed to parse project data from file");
    }

    return QString(); // Success
}

QJsonObject ProjectFile::toJson(const EditableMolecule& molecule,
                               const DrawingStyle& drawingStyle,
                               const PageSettings& pageSettings) {
    QJsonObject json;
    
    // Format version
    json[kFormatVersion] = CurrentFormatVersion;
    
    // Chemistry (molfile)
    QJsonObject chemistry;
    StringResult molResult = molecule.toMolfile();
    if (molResult.success) {
        chemistry[kMolfile] = molResult.value;
    } else {
        chemistry[kMolfile] = QString();
    }
    // Also include name if present
    if (!molecule.name().isEmpty()) {
        chemistry[kName] = molecule.name();
    }
    json[kChemistry] = chemistry;
    
    // Drawing style
    QJsonObject style;
    style[kBondLength] = drawingStyle.bondLength;
    style[kLineWidth] = drawingStyle.lineWidth;
    style[kFontFamily] = drawingStyle.fontFamily;
    style[kFontSize] = drawingStyle.fontSize;
    style[kMonochrome] = drawingStyle.monochrome;
    json[kDrawingStyle] = style;
    
    // Page settings
    QJsonObject page;
    page[kWidthMm] = pageSettings.widthMm;
    page[kHeightMm] = pageSettings.heightMm;
    page[kMarginLeftMm] = pageSettings.marginLeftMm;
    page[kMarginRightMm] = pageSettings.marginRightMm;
    page[kMarginTopMm] = pageSettings.marginTopMm;
    page[kMarginBottomMm] = pageSettings.marginBottomMm;
    json[kPageSettings] = page;
    
    // Fixed physical scale
    json[kModelToMmScale] = drawingStyle.modelToMmScale;
    
    // E/Z intent
    QJsonObject ezIntent;
    for (BondId bondId : molecule.bondIds()) {
        int stereo = molecule.bondStereoDirection(bondId); // Get from molecule
        if (stereo != 0) { // Only save if explicitly set
            ezIntent[QString::number(bondId)] = stereo;
        }
    }
    if (!ezIntent.isEmpty()) {
        json[kEzIntent] = ezIntent;
    }
    
    // S-group state
    QJsonArray sgroups;
    for (SGroupId sgId : molecule.sgroupIds()) {
        QJsonObject sg;
        sg[kId] = sgId;
        sg[kExpanded] = molecule.sgroupExpanded(sgId);
        sg[kLabel] = molecule.sgroupLabel(sgId);
        
        // Get member atoms
        QJsonArray atoms;
        for (AtomId atomId : molecule.sgroupMemberAtomIds(sgId)) {
            atoms.append(atomId);
        }
        sg[kAtoms] = atoms;
        sgroups.append(sg);
    }
    if (!sgroups.isEmpty()) {
        json[kSgroups] = sgroups;
    }
    
    // Annotations (text)
    QJsonArray annotations;
    for (TextId textId : molecule.textAnnotationIds()) {
        double x = 0, y = 0;
        QString content;
        bool bold = false, italic = false;
        if (molecule.textAnnotationContent(textId, x, y, content, bold, italic)) {
            QJsonObject ann;
            ann[kId] = textId;
            ann[kX] = x;
            ann[kY] = y;
            ann[kContent] = content;
            ann[kBold] = bold;
            ann[kItalic] = italic;
            annotations.append(ann);
        }
    }
    if (!annotations.isEmpty()) {
        json[kAnnotations] = annotations;
    }
    
    // Images
    QJsonArray images;
    for (ImageId imageId : molecule.imageIds()) {
        double x = 0, y = 0, w = 0, h = 0;
        QByteArray pngData;
        if (molecule.imageData(imageId, x, y, w, h, pngData)) {
            QJsonObject img;
            img[kId] = imageId;
            img[kX] = x;
            img[kY] = y;
            img[kW] = w;
            img[kH] = h;
            img[kData] = QString::fromLatin1(pngData.toBase64());
            images.append(img);
        }
    }
    if (!images.isEmpty()) {
        json[kImages] = images;
    }
    
    // Reaction arrows
    QJsonArray rxnArrows;
    for (RxnArrowId arrowId : molecule.rxnArrowIds()) {
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        if (molecule.rxnArrowEndpoints(arrowId, x1, y1, x2, y2)) {
            QJsonObject arrow;
            arrow[kId] = arrowId;
            arrow[kX + "1"] = x1;
            arrow[kY + "1"] = y1;
            arrow[kX + "2"] = x2;
            arrow[kY + "2"] = y2;
            arrow[kMode] = molecule.rxnArrowMode(arrowId);
            arrow[kAbove] = molecule.rxnArrowConditionsAbove(arrowId);
            arrow[kBelow] = molecule.rxnArrowConditionsBelow(arrowId);
            
            // Curvature
            double cx = 0, cy = 0;
            if (molecule.rxnArrowCurvature(arrowId, cx, cy)) {
                arrow["curveX"] = cx;
                arrow["curveY"] = cy;
            }
            rxnArrows.append(arrow);
        }
    }
    if (!rxnArrows.isEmpty()) {
        json[kRxnArrows] = rxnArrows;
    }
    
    // Reaction pluses
    QJsonArray rxnPluses;
    for (RxnPlusId plusId : molecule.rxnPlusIds()) {
        double x = 0, y = 0;
        if (molecule.rxnPlusPos(plusId, x, y)) {
            QJsonObject plus;
            plus[kId] = plusId;
            plus[kX] = x;
            plus[kY] = y;
            rxnPluses.append(plus);
        }
    }
    if (!rxnPluses.isEmpty()) {
        json[kRxnPluses] = rxnPluses;
    }
    
    // Multitail arrows
    QJsonArray multitailArrows;
    for (MultitailArrowId mtaId : molecule.multitailArrowIds()) {
        QList<double> points = molecule.multitailArrowPoints(mtaId);
        if (!points.isEmpty()) {
            QJsonObject mta;
            mta[kId] = mtaId;
            QJsonArray pts;
            for (double p : points) {
                pts.append(p);
            }
            mta[kPoints] = pts;
            multitailArrows.append(mta);
        }
    }
    if (!multitailArrows.isEmpty()) {
        json[kMultitailArrows] = multitailArrows;
    }
    
    // R-groups
    QJsonArray rGroups;
    for (int rgroupNumber : molecule.rgroupNumbers()) {
        QString range;
        bool resth = false;
        int ifthen = 0;
        if (molecule.rgroupLogic(rgroupNumber, range, resth, ifthen)) {
            QJsonObject rg;
            rg[kNumber] = rgroupNumber;
            QJsonObject logic;
            logic[kRange] = range;
            logic[kResth] = resth;
            logic[kIfthen] = ifthen;
            rg["logic"] = logic;
            
            QList<int> frags = molecule.rgroupFragmentIds(rgroupNumber);
            QJsonArray fragments;
            for (int frag : frags) {
                fragments.append(frag);
            }
            rg[kFragments] = fragments;
            rGroups.append(rg);
        }
    }
    if (!rGroups.isEmpty()) {
        json[kRGroups] = rGroups;
    }
    
    // Brackets
    QJsonArray brackets;
    int bracketCount = molecule.bracketCount();
    for (int i = 0; i < bracketCount; ++i) {
        double minX = 0, minY = 0, maxX = 0, maxY = 0;
        if (molecule.bracketAt(i, minX, minY, maxX, maxY)) {
            QJsonObject bracket;
            bracket[kMinX] = minX;
            bracket[kMinY] = minY;
            bracket[kMaxX] = maxX;
            bracket[kMaxY] = maxY;
            brackets.append(bracket);
        }
    }
    if (!brackets.isEmpty()) {
        json[kBrackets] = brackets;
    }
    
    // Stereo flags
    QJsonObject stereoFlags;
    stereoFlags["type"] = molecule.stereoFlagsType();
    stereoFlags["groupId"] = molecule.stereoFlagsGroupId();
    json[kStereoFlags] = stereoFlags;
    
    return json;
}

bool ProjectFile::fromJson(const QJsonObject& json,
                          EditableMolecule& molecule,
                          DrawingStyle& drawingStyle,
                          PageSettings& pageSettings) {
    // Validate format version
    if (!json.contains(kFormatVersion)) {
        return false;
    }
    
    int version = json[kFormatVersion].toInt();
    if (version < MinimumSupportedVersion || version > CurrentFormatVersion) {
        return false;
    }
    
    // Load chemistry
    if (json.contains(kChemistry)) {
        QJsonObject chemistry = json[kChemistry].toObject();
        if (chemistry.contains(kMolfile)) {
            QString molfile = chemistry[kMolfile].toString();
            if (!molfile.isEmpty()) {
                // preserveEzIntent=false: this is a genuinely new document replacing
                // whatever was previously loaded into `molecule`, not a same-document
                // write-back round-trip -- per loadFrom's own contract (EditableMolecule.h),
                // passing true here would carry over the PREVIOUSLY open document's
                // per-bond "was E/Z ever specified" flags onto this file's bonds, a stale
                // state leak. E/Z intent for this file is captured fresh from its own
                // molfile stereo encoding, which is what kEzIntent round-trips through.
                //
                // A non-empty molfile that fails to parse means the project
                // file's chemistry is genuinely corrupt -- reporting success
                // with a silently emptied molecule would be worse than
                // rejecting the load outright.
                if (!molecule.loadFrom(molfile, false)) {
                    qWarning() << "ProjectFile::fromJson: Failed to load molfile";
                    return false;
                }
            }
        }
        if (chemistry.contains(kName)) {
            molecule.setName(chemistry[kName].toString());
        }
    }
    
    // Load drawing style
    if (json.contains(kDrawingStyle)) {
        QJsonObject style = json[kDrawingStyle].toObject();
        if (style.contains(kBondLength)) {
            drawingStyle.bondLength = style[kBondLength].toDouble();
        }
        if (style.contains(kLineWidth)) {
            drawingStyle.lineWidth = style[kLineWidth].toDouble();
        }
        if (style.contains(kFontFamily)) {
            drawingStyle.fontFamily = style[kFontFamily].toString();
        }
        if (style.contains(kFontSize)) {
            drawingStyle.fontSize = style[kFontSize].toInt();
        }
        if (style.contains(kMonochrome)) {
            drawingStyle.monochrome = style[kMonochrome].toBool();
        }
        
        // Compute modelToMmScale
        if (json.contains(kModelToMmScale)) {
            drawingStyle.modelToMmScale = json[kModelToMmScale].toDouble();
        } else {
            // Fallback: compute from bondLength
            drawingStyle.modelToMmScale = 5.08 / drawingStyle.bondLength;
        }
    }
    
    // Load page settings
    if (json.contains(kPageSettings)) {
        QJsonObject page = json[kPageSettings].toObject();
        if (page.contains(kWidthMm)) {
            pageSettings.widthMm = page[kWidthMm].toDouble();
        }
        if (page.contains(kHeightMm)) {
            pageSettings.heightMm = page[kHeightMm].toDouble();
        }
        if (page.contains(kMarginLeftMm)) {
            pageSettings.marginLeftMm = page[kMarginLeftMm].toDouble();
        }
        if (page.contains(kMarginRightMm)) {
            pageSettings.marginRightMm = page[kMarginRightMm].toDouble();
        }
        if (page.contains(kMarginTopMm)) {
            pageSettings.marginTopMm = page[kMarginTopMm].toDouble();
        }
        if (page.contains(kMarginBottomMm)) {
            pageSettings.marginBottomMm = page[kMarginBottomMm].toDouble();
        }
    }
    
    // Load E/Z intent
    if (json.contains(kEzIntent)) {
        QJsonObject ezIntent = json[kEzIntent].toObject();
        for (auto it = ezIntent.constBegin(); it != ezIntent.constEnd(); ++it) {
            bool ok;
            BondId bondId = it.key().toInt(&ok);
            if (ok && molecule.findBond(bondId, -1) != -1) { // Check if bond exists
                // Note: E/Z intent is preserved in the molfile itself when loaded
                // So we don't need to explicitly set it here
            }
        }
    }
    
    // Load S-group state
    if (json.contains(kSgroups)) {
        QJsonArray sgroups = json[kSgroups].toArray();
        for (const QJsonValue& val : sgroups) {
            QJsonObject sg = val.toObject();
            SGroupId sgId = sg[kId].toInt();
            if (sg.contains(kExpanded)) {
                molecule.setSGroupExpanded(sgId, sg[kExpanded].toBool());
            }
            if (sg.contains(kLabel)) {
                molecule.setSGroupLabel(sgId, sg[kLabel].toString());
            }
        }
    }
    
    // Load annotations
    if (json.contains(kAnnotations)) {
        QJsonArray annotations = json[kAnnotations].toArray();
        for (const QJsonValue& val : annotations) {
            QJsonObject ann = val.toObject();
            double x = ann[kX].toDouble();
            double y = ann[kY].toDouble();
            QString content = ann[kContent].toString();
            bool bold = ann[kBold].toBool(false);
            bool italic = ann[kItalic].toBool(false);
            int newId = molecule.addTextAnnotation(x, y, content);
            if (newId >= 0 && (bold || italic)) {
                molecule.setTextAnnotation(newId, content, bold, italic);
            }
        }
    }
    
    // Load images
    if (json.contains(kImages)) {
        QJsonArray images = json[kImages].toArray();
        for (const QJsonValue& val : images) {
            QJsonObject img = val.toObject();
            double x = img[kX].toDouble();
            double y = img[kY].toDouble();
            double w = img[kW].toDouble();
            double h = img[kH].toDouble();
            QString dataBase64 = img[kData].toString();
            QByteArray pngData = QByteArray::fromBase64(dataBase64.toLatin1());
            molecule.addImage(x, y, w, h, pngData);
        }
    }
    
    // Load reaction arrows
    if (json.contains(kRxnArrows)) {
        QJsonArray arrows = json[kRxnArrows].toArray();
        for (const QJsonValue& val : arrows) {
            QJsonObject arrow = val.toObject();
            double x1 = arrow[kX + "1"].toDouble();
            double y1 = arrow[kY + "1"].toDouble();
            double x2 = arrow[kX + "2"].toDouble();
            double y2 = arrow[kY + "2"].toDouble();
            RxnArrowId id = molecule.addRxnArrow(x1, y1, x2, y2);
            if (id >= 0) {
                if (arrow.contains(kMode)) {
                    molecule.setRxnArrowMode(id, arrow[kMode].toString());
                }
                if (arrow.contains(kAbove)) {
                    molecule.setRxnArrowConditions(id, arrow[kAbove].toString(), arrow[kBelow].toString());
                }
                if (arrow.contains("curveX") && arrow.contains("curveY")) {
                    molecule.setRxnArrowCurvature(id, arrow["curveX"].toDouble(), arrow["curveY"].toDouble());
                }
            }
        }
    }
    
    // Load reaction pluses
    if (json.contains(kRxnPluses)) {
        QJsonArray pluses = json[kRxnPluses].toArray();
        for (const QJsonValue& val : pluses) {
            QJsonObject plus = val.toObject();
            double x = plus[kX].toDouble();
            double y = plus[kY].toDouble();
            molecule.addRxnPlus(x, y);
        }
    }
    
    // Load multitail arrows
    if (json.contains(kMultitailArrows)) {
        QJsonArray mtas = json[kMultitailArrows].toArray();
        for (const QJsonValue& val : mtas) {
            QJsonObject mta = val.toObject();
            QJsonArray ptsArray = mta[kPoints].toArray();
            QList<double> points;
            for (const QJsonValue& pval : ptsArray) {
                points.append(pval.toDouble());
            }
            if (!points.isEmpty()) {
                molecule.addMultitailArrow(points);
            }
        }
    }
    
    // Load R-groups
    if (json.contains(kRGroups)) {
        QJsonArray rgs = json[kRGroups].toArray();
        for (const QJsonValue& val : rgs) {
            QJsonObject rg = val.toObject();
            int rgroupNumber = rg[kNumber].toInt();
            if (rg.contains("logic")) {
                QJsonObject logic = rg["logic"].toObject();
                QString range = logic[kRange].toString();
                bool resth = logic[kResth].toBool(false);
                int ifthen = logic[kIfthen].toInt(0);
                molecule.setRGroupLogic(rgroupNumber, range, resth, ifthen);
            }
            if (rg.contains(kFragments)) {
                QJsonArray frags = rg[kFragments].toArray();
                for (const QJsonValue& fval : frags) {
                    int fragId = fval.toInt();
                    molecule.addRGroupFragment(rgroupNumber, fragId);
                }
            }
        }
    }
    
    // Load brackets
    if (json.contains(kBrackets)) {
        QJsonArray brs = json[kBrackets].toArray();
        for (const QJsonValue& val : brs) {
            QJsonObject br = val.toObject();
            double minX = br[kMinX].toDouble();
            double minY = br[kMinY].toDouble();
            double maxX = br[kMaxX].toDouble();
            double maxY = br[kMaxY].toDouble();
            molecule.pushBracket(minX, minY, maxX, maxY);
        }
    }
    
    // Load stereo flags
    if (json.contains(kStereoFlags)) {
        QJsonObject stereoFlags = json[kStereoFlags].toObject();
        QString type = stereoFlags["type"].toString();
        int groupId = stereoFlags["groupId"].toInt();
        molecule.setStereoFlagsDocument(type, groupId);
    }
    
    return true;
}

bool ProjectFile::isValidProjectFile(const QString& filePath) {
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return false;
    }
    
    QByteArray data = file.readAll();
    file.close();
    
    QJsonParseError parseError;
    QJsonDocument doc = QJsonDocument::fromJson(data, &parseError);
    
    if (parseError.error != QJsonParseError::NoError) {
        return false;
    }
    
    if (!doc.isObject()) {
        return false;
    }
    
    QJsonObject json = doc.object();
    
    // Must have formatVersion
    if (!json.contains(kFormatVersion)) {
        return false;
    }
    
    int version = json[kFormatVersion].toInt();
    if (version < MinimumSupportedVersion || version > CurrentFormatVersion) {
        return false;
    }
    
    return true;
}

bool ProjectFile::isLegacyRecoveryFile(const QString& filePath) {
    // Check if it's a .ket file
    if (filePath.endsWith(".ket", Qt::CaseInsensitive)) {
        QFile file(filePath);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            return false;
        }
        
        QByteArray data = file.readAll();
        file.close();
        
        // Legacy .ket files are just KET JSON molfiles
        // Check if it looks like a molfile or KET JSON
        QString content = QString::fromUtf8(data);
        return looksLikeMolfile(content) || content.contains("M  END");
    }
    return false;
}

int ProjectFile::getFormatVersion(const QString& filePath) {
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return -1;
    }
    
    QByteArray data = file.readAll();
    file.close();
    
    QJsonParseError parseError;
    QJsonDocument doc = QJsonDocument::fromJson(data, &parseError);
    
    if (parseError.error != QJsonParseError::NoError) {
        return -1;
    }
    
    if (!doc.isObject()) {
        return -1;
    }
    
    QJsonObject json = doc.object();
    
    if (json.contains(kFormatVersion)) {
        return json[kFormatVersion].toInt();
    }
    
    return -1;
}

bool ProjectFile::validateProjectJson(const QJsonObject& json, QString& errorMessage) {
    if (!json.contains(kFormatVersion)) {
        errorMessage = QStringLiteral("Missing formatVersion field");
        return false;
    }
    
    int version = json[kFormatVersion].toInt();
    if (version < MinimumSupportedVersion || version > CurrentFormatVersion) {
        errorMessage = QStringLiteral("Unsupported format version: ") + QString::number(version);
        return false;
    }
    
    // Check for required chemistry field
    if (!json.contains(kChemistry)) {
        errorMessage = QStringLiteral("Missing chemistry field");
        return false;
    }
    
    return true;
}
