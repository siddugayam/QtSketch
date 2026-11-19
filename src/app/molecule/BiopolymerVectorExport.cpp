// src/app/molecule/BiopolymerVectorExport.cpp
#include "BiopolymerVectorExport.h"
#include "BiopolymerState.h"
#include <QDebug>

// ========== Format Conversions ==========

QString vectorExportFormatToString(VectorExportFormat format) {
    switch (format) {
    case VectorExportFormat::Svg: return "SVG";
    case VectorExportFormat::Pdf: return "PDF";
    case VectorExportFormat::Png: return "PNG";
    case VectorExportFormat::Eps: return "EPS";
    case VectorExportFormat::Ps: return "PS";
    default: return "Unknown";
    }
}

VectorExportFormat vectorExportFormatFromString(const QString& formatString) {
    const QString lower = formatString.toLower();
    if (lower == "svg" || lower == "scalable vector graphics") return VectorExportFormat::Svg;
    if (lower == "pdf" || lower == "portable document format") return VectorExportFormat::Pdf;
    if (lower == "png" || lower == "portable network graphics") return VectorExportFormat::Png;
    if (lower == "eps" || lower == "encapsulated postscript") return VectorExportFormat::Eps;
    if (lower == "ps" || lower == "postscript") return VectorExportFormat::Ps;
    return VectorExportFormat::Svg; // Default to SVG
}

// ========== VectorExportOptions ==========

QMap<QString, QVariant> VectorExportOptions::toJson() const {
    QMap<QString, QVariant> json;
    json["format"] = vectorExportFormatToString(format);
    json["includeNumbering"] = includeNumbering;
    json["includeDirection"] = includeDirection;
    json["includeLegends"] = includeLegends;
    json["includeLabels"] = includeLabels;
    json["monochrome"] = monochrome;
    json["colorScheme"] = colorScheme;
    json["width"] = width;
    json["height"] = height;
    json["autoSize"] = autoSize;
    json["margin"] = margin;
    json["fontFamily"] = fontFamily;
    json["fontSize"] = fontSize;
    json["boldLabels"] = boldLabels;
    json["colorByChain"] = colorByChain;
    json["colorByType"] = colorByType;
    json["colorByModification"] = colorByModification;
    return json;
}

VectorExportOptions VectorExportOptions::fromJson(const QMap<QString, QVariant>& json) {
    VectorExportOptions options;
    if (json.contains("format")) {
        options.format = vectorExportFormatFromString(json["format"].toString());
    }
    if (json.contains("includeNumbering")) {
        options.includeNumbering = json["includeNumbering"].toBool();
    }
    if (json.contains("includeDirection")) {
        options.includeDirection = json["includeDirection"].toBool();
    }
    if (json.contains("includeLegends")) {
        options.includeLegends = json["includeLegends"].toBool();
    }
    if (json.contains("includeLabels")) {
        options.includeLabels = json["includeLabels"].toBool();
    }
    if (json.contains("monochrome")) {
        options.monochrome = json["monochrome"].toBool();
    }
    if (json.contains("colorScheme")) {
        options.colorScheme = json["colorScheme"].toString();
    }
    if (json.contains("width")) {
        options.width = json["width"].toInt();
    }
    if (json.contains("height")) {
        options.height = json["height"].toInt();
    }
    if (json.contains("autoSize")) {
        options.autoSize = json["autoSize"].toBool();
    }
    if (json.contains("margin")) {
        options.margin = json["margin"].toInt();
    }
    if (json.contains("fontFamily")) {
        options.fontFamily = json["fontFamily"].toString();
    }
    if (json.contains("fontSize")) {
        options.fontSize = json["fontSize"].toInt();
    }
    if (json.contains("boldLabels")) {
        options.boldLabels = json["boldLabels"].toBool();
    }
    if (json.contains("colorByChain")) {
        options.colorByChain = json["colorByChain"].toBool();
    }
    if (json.contains("colorByType")) {
        options.colorByType = json["colorByType"].toBool();
    }
    if (json.contains("colorByModification")) {
        options.colorByModification = json["colorByModification"].toBool();
    }
    return options;
}

// ========== VectorExportResult ==========

bool VectorExportResult::isSuccess() const {
    return success;
}

bool VectorExportResult::isFailure() const {
    return !success;
}

QByteArray VectorExportResult::toByteArray() const {
    return outputData.toUtf8();
}

// ========== BiopolymerVectorExport ==========

BiopolymerVectorExport::BiopolymerVectorExport(const BiopolymerState& state)
    : m_state(state) {
}

void BiopolymerVectorExport::initialize(int documentId, quint64 documentRevision, 
                                       const QString& documentName) {
    m_documentId = documentId;
    m_documentRevision = documentRevision;
    m_documentName = documentName;
}

int BiopolymerVectorExport::documentId() const {
    return m_documentId;
}

quint64 BiopolymerVectorExport::documentRevision() const {
    return m_documentRevision;
}

QString BiopolymerVectorExport::documentName() const {
    return m_documentName;
}

const BiopolymerState& BiopolymerVectorExport::state() const {
    return m_state;
}

int BiopolymerVectorExport::chainCount() const {
    return m_state.chainCount();
}

int BiopolymerVectorExport::totalResidueCount() const {
    return m_state.totalResidueCount();
}

QList<PolymerType> BiopolymerVectorExport::polymerTypes() const {
    return m_state.polymerTypes();
}

VectorExportResult BiopolymerVectorExport::generateExport(
    VectorExportFormat format, const VectorExportOptions& options) const {
    
    VectorExportResult result;
    result.format = format;
    
    switch (format) {
    case VectorExportFormat::Svg:
        result = generateSvg(options);
        break;
    case VectorExportFormat::Pdf:
        result = generatePdf(options);
        break;
    case VectorExportFormat::Png:
        result = generatePng(options);
        break;
    case VectorExportFormat::Eps:
        result = generateEps(options);
        break;
    case VectorExportFormat::Ps:
    default:
        // Fall back to SVG for unsupported formats
        result = generateSvg(options);
        break;
    }
    
    m_lastResult = result;
    return result;
}

VectorExportResult BiopolymerVectorExport::generateSvg(const VectorExportOptions& options) const {
    VectorExportResult result;
    result.format = VectorExportFormat::Svg;
    
    try {
        // Calculate layouts
        QMap<int, VectorResidue> residues = calculateResidueLayout(options);
        QList<VectorConnection> connections = calculateConnectionLayout(residues, options);
        QList<VectorLegend> legends = calculateLegends(residues, connections, options);
        
        // Calculate total bounds
        QRectF bounds = calculateContentBounds(residues, connections, legends);
        
        // Generate SVG
        QString svg = generateSvgXml(residues, connections, legends, options);
        
        result.success = true;
        result.outputData = svg;
        result.residueCount = residues.size();
        result.connectionCount = connections.size();
        result.chainCount = chainCount();
        result.chainIds = m_state.chainIds();
        result.width = static_cast<int>(bounds.width() + options.margin * 2);
        result.height = static_cast<int>(bounds.height() + options.margin * 2);
        result.contentBounds = bounds;
        
    } catch (const std::exception& e) {
        result.success = false;
        result.errorMessage = QString("SVG generation error: %1").arg(e.what());
    } catch (...) {
        result.success = false;
        result.errorMessage = "Unknown error during SVG generation";
    }
    
    return result;
}

VectorExportResult BiopolymerVectorExport::generatePdf(const VectorExportOptions& options) const {
    VectorExportResult result;
    result.format = VectorExportFormat::Pdf;
    
    // For now, fall back to SVG (PDF generation would require external libraries)
    result = generateSvg(options);
    result.format = VectorExportFormat::Pdf;
    
    return result;
}

VectorExportResult BiopolymerVectorExport::generatePng(const VectorExportOptions& options) const {
    VectorExportResult result;
    result.format = VectorExportFormat::Png;
    
    // PNG generation would require QImage or similar
    // For now, just return SVG as a fallback
    result = generateSvg(options);
    result.format = VectorExportFormat::Png;
    
    return result;
}

VectorExportResult BiopolymerVectorExport::generateEps(const VectorExportOptions& options) const {
    VectorExportResult result;
    result.format = VectorExportFormat::Eps;
    
    // EPS generation would require specialized code
    // For now, fall back to SVG
    result = generateSvg(options);
    result.format = VectorExportFormat::Eps;
    
    return result;
}

QMap<int, VectorResidue> BiopolymerVectorExport::calculateResidueLayout(
    const VectorExportOptions& options) const {
    
    QMap<int, VectorResidue> residues;
    
    // Cell dimensions for layout
    const int cellWidth = 60;
    const int cellHeight = 40;
    const int boxSize = 36;
    const int spacing = 10;
    
    // Track chain positions
    QMap<int, int> chainRow;  // chainId -> row
    QMap<int, int> chainCol;  // chainId -> starting column
    
    // Layout chains horizontally
    int currentRow = 0;
    int currentCol = 0;
    
    for (int chainId : m_state.chainIds()) {
        chainRow[chainId] = currentRow;
        chainCol[chainId] = currentCol;
        currentRow++;
        currentCol = 0;
    }
    
    // Layout residues in each chain
    for (int chainId : m_state.chainIds()) {
        const BiopolymerChain& chain = m_state.chain(chainId);
        int col = chainCol[chainId];
        const int row = chainRow[chainId];
        
        for (int residueId : chain.residueIds()) {
            const BiopolymerResidue& residue = m_state.residue(residueId);
            
            VectorResidue vecRes;
            vecRes.id = residueId;
            vecRes.chainId = chainId;
            vecRes.position = m_state.residuePositionInChain(chainId, residueId);
            vecRes.label = residue.label();
            vecRes.monomerClass = residue.monomerClass();
            
            // Calculate position
            vecRes.center = QPointF(
                col * (cellWidth + spacing) + boxSize / 2,
                row * (cellHeight + spacing) + boxSize / 2
            );
            vecRes.bounds = QRectF(
                col * (cellWidth + spacing),
                row * (cellHeight + spacing),
                boxSize, boxSize
            );
            
            // Get colors
            vecRes.fillColor = getResidueColor(residue, chainId, options);
            vecRes.strokeColor = Qt::black;
            vecRes.strokeWidth = 1;
            
            // Check modifications
            vecRes.isModified = residue.hasModifications();
            vecRes.modifications = residue.modificationNames();
            
            // Check terminals
            vecRes.isNTerminal = m_state.isNTerminal(chainId, residueId);
            vecRes.isCTerminal = m_state.isCTerminal(chainId, residueId);
            vecRes.terminalGroup = residue.terminalGroup();
            
            residues[residueId] = vecRes;
            col++;
        }
    }
    
    return residues;
}

QList<VectorConnection> BiopolymerVectorExport::calculateConnectionLayout(
    const QMap<int, VectorResidue>& residues, 
    const VectorExportOptions& options) const {
    
    QList<VectorConnection> connections;
    
    for (const BiopolymerConnection& conn : m_state.connections()) {
        VectorConnection vecConn;
        vecConn.id = conn.id();
        vecConn.fromId = conn.fromResidueId();
        vecConn.toId = conn.toResidueId();
        vecConn.type = conn.typeString();
        
        // Get positions from residues
        if (residues.contains(conn.fromResidueId()) && residues.contains(conn.toResidueId())) {
            const VectorResidue& fromRes = residues.value(conn.fromResidueId());
            const VectorResidue& toRes = residues.value(conn.toResidueId());
            
            // Simple straight line for now
            vecConn.path.append(fromRes.center);
            vecConn.path.append(toRes.center);
        }
        
        vecConn.strokeColor = getConnectionColor(conn, options);
        vecConn.strokeWidth = 2;
        vecConn.lineStyle = "solid";
        
        // Special styling for certain connection types
        if (conn.isDisulfide()) {
            vecConn.lineStyle = "wavy";
            vecConn.strokeWidth = 3;
        } else if (conn.isInterChain()) {
            vecConn.isDashed = true;
        }
        
        connections.append(vecConn);
    }
    
    return connections;
}

QList<VectorLegend> BiopolymerVectorExport::calculateLegends(
    const QMap<int, VectorResidue>& residues,
    const QList<VectorConnection>& connections,
    const VectorExportOptions& options) const {
    
    QList<VectorLegend> legends;
    
    if (!options.includeLegends) {
        return legends;
    }
    
    // Create chain legend
    if (options.colorByChain && chainCount() > 1) {
        VectorLegend chainLegend;
        chainLegend.title = "Chains";
        
        for (int chainId : m_state.chainIds()) {
            const QString color = getChainColor(chainId, options);
            chainLegend.entries.append(qMakePair(
                QString("Chain %1").arg(chainId),
                color
            ));
        }
        
        legends.append(chainLegend);
    }
    
    // Create modification legend
    bool hasModifications = false;
    for (const VectorResidue& res : residues) {
        if (!res.modifications.isEmpty()) {
            hasModifications = true;
            break;
        }
    }
    
    if (hasModifications) {
        VectorLegend modLegend;
        modLegend.title = "Modifications";
        
        QSet<QString> uniqueMods;
        for (const VectorResidue& res : residues) {
            for (const QString& mod : res.modifications) {
                uniqueMods.insert(mod);
            }
        }
        
        for (const QString& mod : uniqueMods) {
            modLegend.entries.append(qMakePair(mod, "#FF6600")); // Orange for modifications
        }
        
        legends.append(modLegend);
    }
    
    return legends;
}

QRectF BiopolymerVectorExport::calculateContentBounds(
    const QMap<int, VectorResidue>& residues,
    const QList<VectorConnection>& connections,
    const QList<VectorLegend>& legends) const {
    
    QRectF bounds;
    
    // Include residues
    for (const VectorResidue& res : residues) {
        bounds = bounds.united(res.bounds);
    }
    
    // Include connections
    for (const VectorConnection& conn : connections) {
        for (const QPointF& pt : conn.path) {
            bounds = bounds.united(QRectF(pt.x(), pt.y(), 0, 0));
        }
    }
    
    // Add some padding
    bounds.adjust(-10, -10, 10, 10);
    
    return bounds;
}

QString BiopolymerVectorExport::generateSvgXml(
    const QMap<int, VectorResidue>& residues,
    const QList<VectorConnection>& connections,
    const QList<VectorLegend>& legends,
    const VectorExportOptions& options) const {
    
    QStringList svgParts;
    
    // SVG header
    const QRectF bounds = calculateContentBounds(residues, connections, legends);
    const int width = options.width > 0 ? options.width : static_cast<int>(bounds.width() + options.margin * 2);
    const int height = options.height > 0 ? options.height : static_cast<int>(bounds.height() + options.margin * 2);
    
    svgParts << QString("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    svgParts << QString("<svg xmlns=\"http://www.w3.org/2000/svg\" "
                        "xmlns:xlink=\"http://www.w3.org/1999/xlink\" "
                        "width=\"%1\" height=\"%2\" "
                        "viewBox=\"0 0 %1 %2\">").arg(width).arg(height);
    
    // Add defs
    svgParts << generateSvgDefs(options);
    
    // Add connections (behind residues)
    for (const VectorConnection& conn : connections) {
        svgParts << generateSvgConnection(conn, options);
    }
    
    // Add residues
    for (const VectorResidue& res : residues) {
        svgParts << generateSvgResidue(res, options);
    }
    
    // Add direction indicators
    if (options.includeDirection) {
        for (int chainId : m_state.chainIds()) {
            const QList<int> chainResidues = m_state.chain(chainId).residueIds();
            QList<VectorResidue> vecResidues;
            for (int resId : chainResidues) {
                if (residues.contains(resId)) {
                    vecResidues.append(residues.value(resId));
                }
            }
            svgParts << generateDirectionIndicator(chainId, vecResidues, options);
        }
    }
    
    // Add legends
    for (const VectorLegend& legend : legends) {
        svgParts << generateSvgLegend(legend, options);
    }
    
    // Close SVG
    svgParts << "</svg>";
    
    return svgParts.join("\n");
}

QString BiopolymerVectorExport::generateSvgResidue(const VectorResidue& residue, 
                                                  const VectorExportOptions& options) const {
    QStringList parts;
    
    // Use monochrome if requested
    const QString fillColor = options.monochrome ? "#000000" : residue.fillColor;
    const QString strokeColor = options.monochrome ? "#000000" : residue.strokeColor;
    
    // Residue box
    parts << QString("<rect x=\"%1\" y=\"%2\" width=\"%3\" height=\"%4\" "
                     "fill=\"%5\" stroke=\"%6\" stroke-width=\"%7\" "
                     "rx=\"6\" ry=\"6\"/>"
                   ).arg(residue.bounds.x())
                     .arg(residue.bounds.y())
                     .arg(residue.bounds.width())
                     .arg(residue.bounds.height())
                     .arg(fillColor)
                     .arg(strokeColor)
                     .arg(residue.strokeWidth);
    
    // Selection highlight for terminals
    if (residue.isNTerminal || residue.isCTerminal) {
        const QString terminalColor = residue.isNTerminal ? "#00AA00" : "#AA0000";
        parts << QString("<rect x=\"%1\" y=\"%2\" width=\"%3\" height=\"%4\" "
                         "fill=\"none\" stroke=\"%5\" stroke-width=\"2\" "
                         "rx=\"4\" ry=\"4\"/>"
                       ).arg(residue.bounds.x() - 2)
                         .arg(residue.bounds.y() - 2)
                         .arg(residue.bounds.width() + 4)
                         .arg(residue.bounds.height() + 4)
                         .arg(terminalColor);
    }
    
    // Label
    if (options.includeLabels) {
        const QString label = escapeXml(residue.label);
        const int fontSize = options.fontSize > 0 ? options.fontSize : 12;
        const QString fontWeight = options.boldLabels ? "bold" : "normal";
        const QString textColor = options.monochrome ? "#000000" : "#FFFFFF";
        
        parts << QString("<text x=\"%1\" y=\"%2\" "
                         "font-family=\"%3\" font-size=\"%4\" font-weight=\"%5\" "
                         "text-anchor=\"middle\" dominant-baseline=\"central\" "
                         "fill=\"%6\">%7</text>"
                       ).arg(residue.center.x())
                         .arg(residue.center.y())
                         .arg(options.fontFamily.isEmpty() ? "Segoe UI" : escapeXml(options.fontFamily))
                         .arg(fontSize)
                         .arg(fontWeight)
                         .arg(textColor)
                         .arg(label);
    }
    
    // Numbering
    if (options.includeNumbering && residue.position > 0) {
        const QString numberText = QString::number(residue.position);
        const int numFontSize = fontSize - 2;
        
        parts << QString("<text x=\"%1\" y=\"%2\" "
                         "font-family=\"%3\" font-size=\"%4\" "
                         "text-anchor=\"start\" dominant-baseline=\"hanging\" "
                         "fill=\"%5\">%6</text>"
                       ).arg(residue.bounds.x() + 2)
                         .arg(residue.bounds.y() + 2)
                         .arg(options.fontFamily.isEmpty() ? "Segoe UI" : escapeXml(options.fontFamily))
                         .arg(numFontSize)
                         .arg(options.monochrome ? "#000000" : "#000000")
                         .arg(numberText);
    }
    
    // Modification indicator
    if (residue.isModified && options.includeLabels) {
        parts << QString("<circle cx=\"%1\" cy=\"%2\" r=\"4\" "
                         "fill=\"#FF6600\" stroke=\"none\"/>"
                       ).arg(residue.center.x() + residue.bounds.width() / 2 - 8)
                         .arg(residue.center.y() - residue.bounds.height() / 2 + 8);
    }
    
    return parts.join("\n");
}

QString BiopolymerVectorExport::generateSvgConnection(const VectorConnection& connection, 
                                                     const VectorExportOptions& options) const {
    if (connection.path.size() < 2) {
        return "";
    }
    
    const QString color = options.monochrome ? "#000000" : connection.strokeColor;
    const int width = connection.strokeWidth;
    
    // Simple line for now
    const QPointF& start = connection.path.first();
    const QPointF& end = connection.path.last();
    
    if (connection.isDashed) {
        return QString("<line x1=\"%1\" y1=\"%2\" x2=\"%3\" y2=\"%4\" "
                        "stroke=\"%5\" stroke-width=\"%6\" "
                        "stroke-dasharray=\"5,5\" />"
                      ).arg(start.x()).arg(start.y())
                        .arg(end.x()).arg(end.y())
                        .arg(color).arg(width);
    } else if (connection.lineStyle == "wavy") {
        // Simplified wavy line as a zigzag
        const int amplitude = 5;
        const int segments = 10;
        const qreal dx = (end.x() - start.x()) / segments;
        const qreal dy = (end.y() - start.y()) / segments;
        
        QStringList pathParts;
        pathParts << QString("M %1 %2").arg(start.x()).arg(start.y());
        
        for (int i = 1; i < segments; ++i) {
            const qreal x = start.x() + dx * i;
            const qreal y = start.y() + dy * i;
            const qreal offset = (i % 2 == 0) ? amplitude : -amplitude;
            pathParts << QString("L %1 %2").arg(x).arg(y + offset);
        }
        pathParts << QString("L %1 %2").arg(end.x()).arg(end.y());
        
        return QString("<path d=\"%1\" stroke=\"%2\" stroke-width=\"%3\" fill=\"none\"/>")
            .arg(pathParts.join(" ")).arg(color).arg(width);
    } else {
        return QString("<line x1=\"%1\" y1=\"%2\" x2=\"%3\" y2=\"%4\" "
                        "stroke=\"%5\" stroke-width=\"%6\" />"
                      ).arg(start.x()).arg(start.y())
                        .arg(end.x()).arg(end.y())
                        .arg(color).arg(width);
    }
}

QString BiopolymerVectorExport::generateSvgLegend(const VectorLegend& legend, 
                                                 const VectorExportOptions& options) const {
    QStringList parts;
    
    const int legendX = 10;
    const int legendY = 10;
    const int padding = 10;
    const int entryHeight = 20;
    const int entrySpacing = 5;
    
    const QString bgColor = options.monochrome ? "#FFFFFF" : (legend.backgroundColor.isEmpty() ? "#FFFFFF" : legend.backgroundColor);
    const QString textColor = options.monochrome ? "#000000" : (legend.textColor.isEmpty() ? "#000000" : legend.textColor);
    const QString borderColor = legend.borderColor.isEmpty() ? "#000000" : legend.borderColor;
    
    // Calculate legend width and height
    int maxLabelWidth = 0;
    for (const auto& entry : legend.entries) {
        maxLabelWidth = qMax(maxLabelWidth, entry.first.length() * 8);
    }
    
    const int legendWidth = maxLabelWidth + 40 + 20; // label + color box + padding
    const int legendHeight = 30 + legend.entries.size() * (entryHeight + entrySpacing);
    
    // Legend background
    parts << QString("<rect x=\"%1\" y=\"%2\" width=\"%3\" height=\"%4\" "
                     "fill=\"%5\" stroke=\"%6\" stroke-width=\"1\" "
                     "rx=\"4\" ry=\"4\"/>"
                   ).arg(legendX).arg(legendY)
                     .arg(legendWidth).arg(legendHeight)
                     .arg(bgColor).arg(borderColor);
    
    // Legend title
    parts << QString("<text x=\"%1\" y=\"%2\" "
                     "font-family=\"Segoe UI\" font-size=\"12\" font-weight=\"bold\" "
                     "text-anchor=\"start\" dominant-baseline=\"hanging\" "
                     "fill=\"%3\">%4</text>"
                   ).arg(legendX + padding)
                     .arg(legendY + padding)
                     .arg(textColor)
                     .arg(escapeXml(legend.title));
    
    // Legend entries
    int currentY = legendY + padding + 20;
    for (const auto& entry : legend.entries) {
        // Color box
        parts << QString("<rect x=\"%1\" y=\"%2\" width=\"16\" height=\"16\" "
                         "fill=\"%3\" stroke=\"none\"/>"
                       ).arg(legendX + padding)
                         .arg(currentY)
                         .arg(entry.second);
        
        // Label
        parts << QString("<text x=\"%1\" y=\"%2\" "
                         "font-family=\"Segoe UI\" font-size=\"10\" "
                         "text-anchor=\"start\" dominant-baseline=\"central\" "
                         "fill=\"%3\">%4</text>"
                       ).arg(legendX + padding + 22)
                         .arg(currentY + 8)
                         .arg(textColor)
                         .arg(escapeXml(entry.first));
        
        currentY += entryHeight + entrySpacing;
    }
    
    return parts.join("\n");
}

QString BiopolymerVectorExport::generateSvgDefs(const VectorExportOptions& options) const {
    QStringList defs;
    
    defs << "<defs>";
    
    // Add filters or gradients if needed
    if (options.monochrome) {
        // Monochrome doesn't need special defs
    }
    
    defs << "</defs>";
    
    return defs.join("\n");
}

QString BiopolymerVectorExport::generateDirectionIndicator(
    int chainId, const QList<VectorResidue>& chainResidues, 
    const VectorExportOptions& options) const {
    
    if (chainResidues.isEmpty()) {
        return "";
    }
    
    const VectorResidue& first = chainResidues.first();
    const VectorResidue& last = chainResidues.last();
    
    // Determine direction based on chain type
    const BiopolymerChain& chain = m_state.chain(chainId);
    const QString directionText = (chain.polymerType() == PolymerType::Protein) ? 
        "N→C" : "5'→3'";
    
    // Position arrow at the start of the chain
    const QPointF arrowStart = QPointF(
        first.bounds.x() - 20,
        first.center.y()
    );
    
    const QPointF arrowEnd = QPointF(
        first.bounds.x() - 5,
        first.center.y()
    );
    
    const QString color = options.monochrome ? "#000000" : "#0066CC";
    
    // Draw arrow line
    QStringList parts;
    parts << QString("<line x1=\"%1\" y1=\"%2\" x2=\"%3\" y2=\"%4\" "
                     "stroke=\"%5\" stroke-width=\"2\" "
                     "marker-end=\"url(#arrowhead)\" />"
                   ).arg(arrowStart.x()).arg(arrowStart.y())
                     .arg(arrowEnd.x()).arg(arrowEnd.y())
                     .arg(color);
    
    // Draw direction text
    parts << QString("<text x=\"%1\" y=\"%2\" "
                     "font-family=\"Segoe UI\" font-size=\"10\" "
                     "text-anchor=\"end\" dominant-baseline=\"central\" "
                     "fill=\"%3\">%4</text>"
                   ).arg(arrowStart.x() - 5)
                     .arg(arrowStart.y() - 5)
                     .arg(color)
                     .arg(directionText);
    
    return parts.join("\n");
}

QString BiopolymerVectorExport::generateSvgPath(const QList<QPointF>& points, bool curved) const {
    if (points.size() < 2) {
        return "";
    }
    
    QStringList pathParts;
    pathParts << QString("M %1 %2").arg(points.first().x()).arg(points.first().y());
    
    if (curved) {
        // Create smooth curve
        for (int i = 1; i < points.size() - 1; ++i) {
            const QPointF& p0 = points[i - 1];
            const QPointF& p1 = points[i];
            const QPointF& p2 = points[i + 1];
            
            const QPointF cp1 = QPointF(
                p1.x() + (p2.x() - p0.x()) / 6,
                p1.y() + (p2.y() - p0.y()) / 6
            );
            
            pathParts << QString("S %1 %2 %3 %4")
                .arg(cp1.x()).arg(cp1.y())
                .arg(p2.x()).arg(p2.y());
        }
    } else {
        // Straight lines
        for (int i = 1; i < points.size(); ++i) {
            pathParts << QString("L %1 %2").arg(points[i].x()).arg(points[i].y());
        }
    }
    
    return pathParts.join(" ");
}

QString BiopolymerVectorExport::generateSvgArrowhead(const QPointF& start, const QPointF& end, 
                                                     const QString& color, int width) const {
    const qreal dx = end.x() - start.x();
    const qreal dy = end.y() - start.y();
    const qreal length = qSqrt(dx * dx + dy * dy);
    
    if (length < 1) {
        return "";
    }
    
    const qreal angle = qAtan2(dy, dx);
    const qreal arrowLength = width * 2;
    const qreal arrowAngle = qPi / 6; // 30 degrees
    
    const QPointF arrow1 = QPointF(
        end.x() - arrowLength * qCos(angle - arrowAngle),
        end.y() - arrowLength * qSin(angle - arrowAngle)
    );
    const QPointF arrow2 = QPointF(
        end.x() - arrowLength * qCos(angle + arrowAngle),
        end.y() - arrowLength * qSin(angle + arrowAngle)
    );
    
    return QString("<path d=\"M %1 %2 L %3 %4 L %5 %6 Z\" "
                    "fill=\"%7\" stroke=\"none\"/>")
        .arg(end.x()).arg(end.y())
        .arg(arrow1.x()).arg(arrow1.y())
        .arg(arrow2.x()).arg(arrow2.y())
        .arg(color);
}

QString BiopolymerVectorExport::getResidueColor(const BiopolymerResidue& residue, 
                                                 int chainId, 
                                                 const VectorExportOptions& options) const {
    // Monochrome mode
    if (options.monochrome) {
        return "#000000";
    }
    
    // Color by chain
    if (options.colorByChain) {
        return getChainColor(chainId, options);
    }
    
    // Color by type
    if (options.colorByType) {
        const PolymerType type = m_state.chain(chainId).polymerType();
        switch (type) {
        case PolymerType::Protein:
            return "#4CAF50"; // Green
        case PolymerType::Dna:
        case PolymerType::Rna:
            return "#2196F3"; // Blue
        default:
            return "#9E9E9E"; // Gray
        }
    }
    
    // Color by modification
    if (options.colorByModification && residue.hasModifications()) {
        return "#FF9800"; // Orange for modified
    }
    
    // Default color based on monomer class
    if (residue.monomerClass() == "Base") {
        return "#8764B8"; // Purple
    }
    
    // Default to accent color
    return "#0067C0";
}

QString BiopolymerVectorExport::getConnectionColor(const BiopolymerConnection& connection, 
                                                    const VectorExportOptions& options) const {
    if (options.monochrome) {
        return "#000000";
    }
    
    if (connection.isDisulfide()) {
        return "#FFC107"; // Amber for disulfide
    }
    
    if (connection.isInterChain()) {
        return "#607D8B"; // Blue-gray for inter-chain
    }
    
    return "#212121"; // Dark gray
}

QString BiopolymerVectorExport::getChainColor(int chainId, 
                                               const VectorExportOptions& options) const {
    // Check for custom chain color
    if (options.chainColors.contains(chainId)) {
        return options.chainColors[chainId];
    }
    
    // Default chain colors
    static const QList<QString> defaultChainColors = {
        "#0067C0", // Blue
        "#E81123", // Red
        "#00AA44", // Green
        "#F58220", // Orange
        "#8764B8", // Purple
        "#00AA88", // Teal
        "#FF0084", // Pink
        "#0088CC", // Light Blue
    };
    
    const int index = chainId % defaultChainColors.size();
    if (index >= 0 && index < defaultChainColors.size()) {
        return defaultChainColors[index];
    }
    
    return "#9E9E9E"; // Gray fallback
}

QString BiopolymerVectorExport::getMonochromeColor(int index, int total) const {
    Q_UNUSED(index);
    Q_UNUSED(total);
    return "#000000";
}

int BiopolymerVectorExport::getResiduePositionInChain(int chainId, int residueId) const {
    return m_state.residuePositionInChain(chainId, residueId);
}

bool BiopolymerVectorExport::isResidueAtNTerminal(int chainId, int residueId) const {
    return m_state.isNTerminal(chainId, residueId);
}

bool BiopolymerVectorExport::isResidueAtCTerminal(int chainId, int residueId) const {
    return m_state.isCTerminal(chainId, residueId);
}

QList<int> BiopolymerVectorExport::getChainResiduesInOrder(int chainId) const {
    return m_state.chain(chainId).residueIds();
}

// ========== Static Helpers ==========

QString BiopolymerVectorExport::escapeXml(const QString& text) {
    QString result = text;
    result.replace("&", "&amp;");
    result.replace("<", "&lt;");
    result.replace(">", "&gt;");
    result.replace("\"", "&quot;");
    result.replace("'", "&apos;");
    return result;
}

QString BiopolymerVectorExport::formatSvgStyle(const QMap<QString, QString>& properties) {
    QStringList styleParts;
    for (auto it = properties.constBegin(); it != properties.constEnd(); ++it) {
        styleParts << QString("%1:%2").arg(it.key()).arg(it.value());
    }
    return styleParts.join(";");
}

QString BiopolymerVectorExport::formatMimeType(VectorExportFormat format) {
    switch (format) {
    case VectorExportFormat::Svg: return "image/svg+xml";
    case VectorExportFormat::Pdf: return "application/pdf";
    case VectorExportFormat::Png: return "image/png";
    case VectorExportFormat::Eps: return "application/postscript";
    case VectorExportFormat::Ps: return "application/postscript";
    default: return "image/svg+xml";
    }
}

QString BiopolymerVectorExport::formatFileExtension(VectorExportFormat format) {
    switch (format) {
    case VectorExportFormat::Svg: return ".svg";
    case VectorExportFormat::Pdf: return ".pdf";
    case VectorExportFormat::Png: return ".png";
    case VectorExportFormat::Eps: return ".eps";
    case VectorExportFormat::Ps: return ".ps";
    default: return ".svg";
    }
}

QString BiopolymerVectorExport::formatDisplayName(VectorExportFormat format) {
    switch (format) {
    case VectorExportFormat::Svg: return "SVG";
    case VectorExportFormat::Pdf: return "PDF";
    case VectorExportFormat::Png: return "PNG";
    case VectorExportFormat::Eps: return "EPS";
    case VectorExportFormat::Ps: return "PostScript";
    default: return "Vector";
    }
}

QString BiopolymerVectorExport::formatDescription(VectorExportFormat format) {
    switch (format) {
    case VectorExportFormat::Svg: return "Scalable Vector Graphics";
    case VectorExportFormat::Pdf: return "Portable Document Format";
    case VectorExportFormat::Png: return "Portable Network Graphics";
    case VectorExportFormat::Eps: return "Encapsulated PostScript";
    case VectorExportFormat::Ps: return "PostScript";
    default: return "Vector Graphics";
    }
}

bool BiopolymerVectorExport::formatIsVectorBased(VectorExportFormat format) {
    return format != VectorExportFormat::Png;
}

VectorExportOptions BiopolymerVectorExport::defaultOptions(VectorExportFormat format) {
    Q_UNUSED(format);
    VectorExportOptions options;
    options.format = format;
    options.includeNumbering = true;
    options.includeDirection = true;
    options.includeLegends = true;
    options.includeLabels = true;
    options.monochrome = false;
    options.width = 800;
    options.height = 600;
    options.autoSize = true;
    options.margin = 20;
    options.fontSize = 12;
    options.boldLabels = false;
    options.colorByChain = true;
    return options;
}

// ========== VectorResidue ==========

bool operator==(const VectorResidue& lhs, const VectorResidue& rhs) {
    return lhs.id == rhs.id && lhs.chainId == rhs.chainId;
}

// ========== VectorConnection ==========

bool operator==(const VectorConnection& lhs, const VectorConnection& rhs) {
    return lhs.id == rhs.id;
}
