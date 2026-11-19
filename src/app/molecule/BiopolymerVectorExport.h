// src/app/molecule/BiopolymerVectorExport.h
#ifndef BIOPOLYMERVECTOREXPORT_H
#define BIOPOLYMERVECTOREXPORT_H

#include "BiopolymerState.h"
#include <QString>
#include <QList>
#include <QMap>
#include <QRectF>
#include <QPointF>

// ========== Vector Export Format ==========

enum class VectorExportFormat {
    Svg,           // Scalable Vector Graphics
    Pdf,           // Portable Document Format
    Png,           // Portable Network Graphics (raster, but useful for quick previews)
    Eps,           // Encapsulated PostScript
    Ps             // PostScript
};

// Convert vector export format to string
QString vectorExportFormatToString(VectorExportFormat format);

// Parse vector export format from string
VectorExportFormat vectorExportFormatFromString(const QString& formatString);

// ========== Vector Export Options ==========

// Options for vector export
struct VectorExportOptions {
    VectorExportFormat format = VectorExportFormat::Svg;
    
    // Layout options
    bool includeNumbering = true;        // Include residue numbers
    bool includeDirection = true;       // Include direction indicators (N->C, 5'->3')
    bool includeLegends = true;          // Include legend boxes
    bool includeLabels = true;          // Include residue labels
    
    // Styling options
    bool monochrome = false;             // Monochrome output (black/white only)
    QString colorScheme;                // Color scheme name ("default", "cpk", "amino", etc.)
    
    // Sizing options
    int width = 800;                    // Output width in pixels (for raster) or points
    int height = 600;                   // Output height in pixels or points
    bool autoSize = true;               // Automatically determine size from content
    int margin = 20;                    // Margin around content in pixels/points
    
    // Font options
    QString fontFamily;                 // Font family for labels
    int fontSize = 12;                  // Base font size
    bool boldLabels = false;            // Use bold for labels
    
    // Chain coloring
    bool colorByChain = true;           // Color residues by chain
    QMap<int, QString> chainColors;      // Custom colors per chain
    
    // Residue coloring
    bool colorByType = false;           // Color by residue type
    bool colorByModification = false;   // Color by modification state
    
    // Serialization
    QMap<QString, QVariant> toJson() const;
    static VectorExportOptions fromJson(const QMap<QString, QVariant>& json);
};

// ========== Vector Residue Element ==========

// Represents a residue in vector output
struct VectorResidue {
    int id = -1;
    int chainId = -1;
    int position = -1;                  // Position in chain (1-indexed)
    QString label;                     // Residue label (e.g., "A", "T", "G", "C", "ALA")
    QString monomerClass;               // Monomer class (Base, AminoAcid, etc.)
    QPointF center;                    // Center position in output coordinates
    QRectF bounds;                     // Bounding box in output coordinates
    
    // Visual properties
    QString fillColor;
    QString strokeColor;
    int strokeWidth = 1;
    
    // Modification state
    bool isModified = false;
    QStringList modifications;         // List of modification names
    
    // Terminal state
    bool isNTerminal = false;
    bool isCTerminal = false;
    QString terminalGroup;             // Terminal group name if applicable
};

// ========== Vector Connection Element ==========

// Represents a connection in vector output
struct VectorConnection {
    int id = -1;
    int fromId = -1;                    // Source residue ID
    int toId = -1;                      // Destination residue ID
    QString type;                      // Connection type (peptide, phosphodiester, disulfide, etc.)
    QList<QPointF> path;                // Path points for drawing the connection
    
    // Visual properties
    QString strokeColor;
    int strokeWidth = 2;
    bool isDashed = false;              // Dashed line for special connections
    QString lineStyle;                 // Line style (solid, dashed, dotted, wavy)
};

// ========== Vector Legend Element ==========

// Represents a legend entry in vector output
struct VectorLegend {
    QString title;
    QList<QPair<QString, QString>> entries;  // Pair of (label, color/description)
    QRectF bounds;                     // Position and size of legend box
    
    // Visual properties
    QString backgroundColor;
    QString textColor;
    QString borderColor;
    int borderWidth = 1;
};

// ========== Vector Export Result ==========

// Result of a vector export operation
struct VectorExportResult {
    VectorExportFormat format = VectorExportFormat::Svg;
    bool success = false;
    QString outputData;                 // The generated vector data (SVG XML, PDF binary, etc.)
    QString errorMessage;
    
    // Metadata
    int width = 0;
    int height = 0;
    QRectF contentBounds;               // Bounds of actual content (excluding margins)
    
    // Content summary
    int residueCount = 0;
    int connectionCount = 0;
    int chainCount = 0;
    QList<int> chainIds;
    
    // Check if export succeeded
    bool isSuccess() const;
    
    // Check if export failed
    bool isFailure() const;
    
    // Get output as bytes (for binary formats)
    QByteArray toByteArray() const;
};

// ========== BiopolymerVectorExport ==========

// Class for exporting biopolymer structures to vector formats
//
// This class generates vector graphics (SVG, PDF, etc.) from biopolymer
// structures with support for:
// - Residue numbering
// - Direction indicators (N->C, 5'->3')
// - Legends
// - Monochrome output
// - Custom styling
class BiopolymerVectorExport {
public:
    // Constructor
    explicit BiopolymerVectorExport(const BiopolymerState& state);
    
    // Initialize with document information
    void initialize(int documentId, quint64 documentRevision, const QString& documentName);
    
    // ========== Export Generation ==========
    
    // Generate vector export for a specific format
    VectorExportResult generateExport(VectorExportFormat format, 
                                     const VectorExportOptions& options = VectorExportOptions()) const;
    
    // Generate SVG export
    VectorExportResult generateSvg(const VectorExportOptions& options = VectorExportOptions()) const;
    
    // Generate PDF export (requires external library or Qt PDF support)
    VectorExportResult generatePdf(const VectorExportOptions& options = VectorExportOptions()) const;
    
    // Generate PNG export (raster, for quick preview)
    VectorExportResult generatePng(const VectorExportOptions& options = VectorExportOptions()) const;
    
    // Generate EPS export
    VectorExportResult generateEps(const VectorExportOptions& options = VectorExportOptions()) const;
    
    // ========== Layout Calculation ==========
    
    // Calculate layout for all residues
    QMap<int, VectorResidue> calculateResidueLayout(const VectorExportOptions& options) const;
    
    // Calculate layout for all connections
    QList<VectorConnection> calculateConnectionLayout(
        const QMap<int, VectorResidue>& residues, 
        const VectorExportOptions& options) const;
    
    // Calculate layout for legends
    QList<VectorLegend> calculateLegends(
        const QMap<int, VectorResidue>& residues,
        const QList<VectorConnection>& connections,
        const VectorExportOptions& options) const;
    
    // Calculate total bounds for content
    QRectF calculateContentBounds(
        const QMap<int, VectorResidue>& residues,
        const QList<VectorConnection>& connections,
        const QList<VectorLegend>& legends) const;
    
    // ========== SVG Generation ==========
    
    // Generate SVG XML for the entire structure
    QString generateSvgXml(const QMap<int, VectorResidue>& residues,
                           const QList<VectorConnection>& connections,
                           const QList<VectorLegend>& legends,
                           const VectorExportOptions& options) const;
    
    // Generate SVG for a single residue
    QString generateSvgResidue(const VectorResidue& residue, const VectorExportOptions& options) const;
    
    // Generate SVG for a single connection
    QString generateSvgConnection(const VectorConnection& connection, const VectorExportOptions& options) const;
    
    // Generate SVG for a legend
    QString generateSvgLegend(const VectorLegend& legend, const VectorExportOptions& options) const;
    
    // Generate SVG defs (styles, gradients, etc.)
    QString generateSvgDefs(const VectorExportOptions& options) const;
    
    // ========== Direction Indicators ==========
    
    // Generate direction indicator for a chain
    QString generateDirectionIndicator(int chainId, 
                                       const QList<VectorResidue>& chainResidues,
                                       const VectorExportOptions& options) const;
    
    // ========== Numbering ==========
    
    // Generate numbering for residues
    QString generateNumbering(const VectorResidue& residue, const VectorExportOptions& options) const;
    
    // ========== Color Management ==========
    
    // Get color for a residue based on options
    QString getResidueColor(const BiopolymerResidue& residue, int chainId, 
                            const VectorExportOptions& options) const;
    
    // Get color for a connection based on options
    QString getConnectionColor(const BiopolymerConnection& connection, 
                               const VectorExportOptions& options) const;
    
    // Get color for a chain
    QString getChainColor(int chainId, const VectorExportOptions& options) const;
    
    // Get monochrome color for a feature
    QString getMonochromeColor(int index, int total) const;
    
    // ========== Utility Methods ==========
    
    // Get document information
    int documentId() const;
    quint64 documentRevision() const;
    QString documentName() const;
    
    // Get biopolymer state
    const BiopolymerState& state() const;
    
    // Get summary information
    int chainCount() const;
    int totalResidueCount() const;
    QList<PolymerType> polymerTypes() const;
    
    // ========== Static Helpers ==========
    
    // Escape XML special characters
    static QString escapeXml(const QString& text);
    
    // Format SVG style string
    static QString formatSvgStyle(const QMap<QString, QString>& properties);
    
    // Get format MIME type
    static QString formatMimeType(VectorExportFormat format);
    
    // Get format file extension
    static QString formatFileExtension(VectorExportFormat format);
    
    // Get format display name
    static QString formatDisplayName(VectorExportFormat format);
    
    // Get format description
    static QString formatDescription(VectorExportFormat format);
    
    // Check if format is vector-based
    static bool formatIsVectorBased(VectorExportFormat format);
    
    // Get default options for format
    static VectorExportOptions defaultOptions(VectorExportFormat format);

private:
    // ========== Private Members ==========
    
    const BiopolymerState& m_state;
    
    int m_documentId = -1;
    quint64 m_documentRevision = 0;
    QString m_documentName;
    
    mutable VectorExportResult m_lastResult;
    
    // ========== Private Helper Methods ==========
    
    // Get residue position in chain
    int getResiduePositionInChain(int chainId, int residueId) const;
    
    // Check if residue is at chain terminal
    bool isResidueAtNTerminal(int chainId, int residueId) const;
    bool isResidueAtCTerminal(int chainId, int residueId) const;
    
    // Get chain residues in order
    QList<int> getChainResiduesInOrder(int chainId) const;
    
    // Calculate residue position in layout
    QPointF calculateResiduePosition(int chainId, int residueId, 
                                     const VectorExportOptions& options) const;
    
    // Generate SVG path for a connection
    QString generateSvgPath(const QList<QPointF>& points, bool curved = false) const;
    
    // Generate SVG arrowhead
    QString generateSvgArrowhead(const QPointF& start, const QPointF& end, 
                                const QString& color, int width = 10) const;
};

#endif // BIOPOLYMERVECTOREXPORT_H
