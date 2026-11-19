// src/app/molecule/IndigoMonomerAdapter.h
#ifndef INDIGOMONOMERADAPTER_H
#define INDIGOMONOMERADAPTER_H

// Monomer adapter over the bundled Indigo library.
//
// Provides concrete monomer definitions, stereochemistry, attachments, and terminal chemistry
// for standard 20 peptide residues and canonical DNA/RNA bases.
//
// Produces explicit expansion correspondence between residue instances and generated atoms.
// Validates correspondence before accepting generated chemistry.

#include <QString>
#include <QHash>
#include <QList>
#include <QMap>
#include <QSet>

using AtomId = int;
using BondId = int;

// Forward declaration
class EditableMolecule;

// Polymer type (matches BiopolymerState)
enum class PolymerType {
    Peptide,
    DNA,
    RNA,
    Unknown
};

// Monomer template definition
struct MonomerTemplate {
    QString symbol;               // Single-letter code
    QString alias;                // Full name
    QString monomerClass;         // "AminoAcid" or "Base"
    bool ambiguous = false;
    
    // Indigo SMILES for this monomer (without attachment points)
    QString smiles;
    
    // Attachment points for polymer connections
    // For peptides: 0 = N-terminus, 1 = C-terminus
    // For nucleotides: 0 = 5' phosphate, 1 = 3' hydroxyl, 2 = base
    QList<int> attachmentPoints;
    
    // Atomic details for atom mapping
    QList<QString> atomLabels;    // Expected atom labels in order
    QList<int> atomCharges;
    QList<int> atomValencies;
    
    // For peptides: stereochemistry at chiral centers
    QString stereochemistry;      // "L" or "D"
    QList<int> chiralAtomIndices; // Indices of chiral atoms (0-based in atomLabels)
    
    // For nucleotides: base, sugar, phosphate components
    QList<int> baseAtomIndices;
    QList<int> sugarAtomIndices;
    QList<int> phosphateAtomIndices;
};

// Terminal group template
struct TerminalGroupTemplate {
    QString name;
    QString smiles;
    QList<int> attachmentIndices;  // Where this attaches to residue
};

// Expansion result: maps a residue to generated atoms
struct ExpansionResult {
    bool success = false;
    QString errorMessage;
    
    // Generated structure info
    QString molfile;              // MOL file of the expanded monomer
    QList<AtomId> createdAtomIds; // IDs of created atoms (from EditableMolecule)
    QList<BondId> createdBondIds; // IDs of created bonds (from EditableMolecule)
    
    // Atom mapping: position in template -> AtomId
    QMap<int, AtomId> templateIndexToAtomId;
    
    // Attachment point AtomIds
    QMap<int, AtomId> attachmentPointAtomIds;
};

// Connection result: maps a bond between residues
struct ConnectionResult {
    bool success = false;
    QString errorMessage;
    
    BondId createdBondId = -1;
    AtomId fromAtomId = -1;
    AtomId toAtomId = -1;
};

// The main adapter class
class IndigoMonomerAdapter {
public:
    IndigoMonomerAdapter();
    ~IndigoMonomerAdapter();

    // ========== Monomer Template Access ==========
    
    // Get template for a monomer symbol and polymer type
    const MonomerTemplate* templateForSymbol(const QString& symbol, PolymerType type) const;
    
    // Check if a symbol is valid for a polymer type
    bool isValidSymbol(const QString& symbol, PolymerType type) const;
    
    // Get all valid symbols for a polymer type
    QList<QString> validSymbols(PolymerType type) const;
    
    // Get all templates for a polymer type
    QList<MonomerTemplate> templates(PolymerType type) const;
    
    // ========== Monomer Expansion ==========
    
    // Expand a monomer to atoms using EditableMolecule
    // Returns ExpansionResult with mapping information
    ExpansionResult expandMonomer(EditableMolecule& molecule, 
                                   const QString& symbol, 
                                   PolymerType type,
                                   double x, double y) const;
    
    // Expand a monomer with specific atom labels
    ExpansionResult expandMonomerWithLabels(EditableMolecule& molecule,
                                              const QString& symbol,
                                              PolymerType type,
                                              const QMap<int, QString>& atomLabelOverrides,
                                              double x, double y) const;
    
    // ========== Terminal Group Expansion ==========
    
    // Expand N-terminus (peptide)
    ExpansionResult expandNTerminus(EditableMolecule& molecule, AtomId attachTo, double x, double y) const;
    
    // Expand C-terminus (peptide)
    ExpansionResult expandCTerminus(EditableMolecule& molecule, AtomId attachTo, double x, double y) const;
    
    // Expand 5' phosphate (nucleotide)
    ExpansionResult expandPhosphate5Prime(EditableMolecule& molecule, AtomId attachTo, double x, double y) const;
    
    // Expand 3' phosphate (nucleotide)
    ExpansionResult expandPhosphate3Prime(EditableMolecule& molecule, AtomId attachTo, double x, double y) const;
    
    // ========== Residue Connection ==========
    
    // Connect two residues with a bond
    ConnectionResult connectResidues(EditableMolecule& molecule,
                                      const QString& fromSymbol, PolymerType fromType,
                                      const QString& toSymbol, PolymerType toType,
                                      AtomId fromAtomId, AtomId toAtomId) const;
    
    // Connect using attachment point indices
    ConnectionResult connectByAttachmentPoints(EditableMolecule& molecule,
                                               const QString& fromSymbol, PolymerType fromType,
                                               const QString& toSymbol, PolymerType toType,
                                               int fromAttachmentIndex,
                                               int toAttachmentIndex,
                                               AtomId fromAtomId, AtomId toAtomId) const;
    
    // ========== Validation ==========
    
    // Validate that expanded atoms match the template
    bool validateExpansion(const EditableMolecule& molecule,
                            const ExpansionResult& result,
                            const QString& symbol,
                            PolymerType type) const;
    
    // Validate a bond connection between residues
    bool validateConnection(const EditableMolecule& molecule,
                             const ConnectionResult& result,
                             const QString& fromSymbol, PolymerType fromType,
                             const QString& toSymbol, PolymerType toType) const;
    
    // ========== Utility Methods ==========
    
    // Get the attachment point atom for a residue
    AtomId getAttachmentPointAtom(const EditableMolecule& molecule,
                                  const QString& symbol,
                                  PolymerType type,
                                  int attachmentIndex,
                                  const QMap<int, AtomId>& templateIndexToAtomId) const;
    
    // Check if a residue has a specific modification
    bool hasModification(const QString& symbol, PolymerType type, const QString& modification) const;
    
    // Get stereochemistry info for a residue
    QString getStereochemistry(const QString& symbol, PolymerType type) const;
    
    // Get the expected number of heavy atoms for a residue
    int expectedHeavyAtomCount(const QString& symbol, PolymerType type) const;
    
    // Get atom label at a specific position in the template
    QString atomLabelAtPosition(const QString& symbol, PolymerType type, int position) const;
    
    // ========== Static Helper Methods ==========
    
    // Get string representation of polymer type
    static QString polymerTypeToString(PolymerType type);
    
    // Normalize polymer type string
    static QString normalizePolymerTypeString(const QString& type);
    
    // Normalize monomer symbol
    static QString normalizeMonomerSymbol(const QString& symbol);
    
    // Check if a polymer type supports a modification
    static bool supportsModification(PolymerType type, const QString& modification);

private:
    // ========== Template Data ==========
    
    // Standard 20 amino acids (L-form by default)
    static const QList<MonomerTemplate> s_peptideTemplates;
    
    // DNA bases
    static const QList<MonomerTemplate> s_dnaTemplates;
    
    // RNA bases
    static const QList<MonomerTemplate> s_rnaTemplates;
    
    // Ambiguous amino acids
    static const QList<MonomerTemplate> s_ambiguousPeptideTemplates;
    
    // Ambiguous nucleotides
    static const QList<MonomerTemplate> s_ambiguousNucleotideTemplates;
    
    // Terminal group templates
    static const QList<TerminalGroupTemplate> s_peptideTerminalGroups;
    static const QList<TerminalGroupTemplate> s_nucleotideTerminalGroups;
    
    // ========== Private Methods ==========
    
    // Create peptide templates
    static QList<MonomerTemplate> createPeptideTemplates();
    
    // Create DNA templates
    static QList<MonomerTemplate> createDnaTemplates();
    
    // Create RNA templates
    static QList<MonomerTemplate> createRnaTemplates();
    
    // Create ambiguous peptide templates
    static QList<MonomerTemplate> createAmbiguousPeptideTemplates();
    
    // Create ambiguous nucleotide templates
    static QList<MonomerTemplate> createAmbiguousNucleotideTemplates();
    
    // Create terminal group templates
    static QList<TerminalGroupTemplate> createPeptideTerminalGroups();
    static QList<TerminalGroupTemplate> createNucleotideTerminalGroups();
    
    // Build full template map (including ambiguous)
    static QHash<QString, MonomerTemplate> buildFullTemplateMap(PolymerType type);
    
    // Get templates for a specific polymer type
    const QList<MonomerTemplate>& templatesForType(PolymerType type) const;
    
    // Get ambiguous templates for a specific polymer type
    const QList<MonomerTemplate>& ambiguousTemplatesForType(PolymerType type) const;
    
    // Parse Indigo molfile to extract atom/bond IDs
    QMap<int, AtomId> parseMolfileAtomMap(const QString& molfile, 
                                         const EditableMolecule& molecule,
                                         const QMap<int, AtomId>& existingMap) const;
    
    // Get the template map for a polymer type
    QHash<QString, MonomerTemplate> getTemplateMap(PolymerType type) const;
    
    // Member variables
    QHash<QString, MonomerTemplate> m_peptideTemplateMap;
    QHash<QString, MonomerTemplate> m_dnaTemplateMap;
    QHash<QString, MonomerTemplate> m_rnaTemplateMap;
    
    QHash<QString, MonomerTemplate> m_fullPeptideTemplateMap;
    QHash<QString, MonomerTemplate> m_fullDnaTemplateMap;
    QHash<QString, MonomerTemplate> m_fullRnaTemplateMap;
};

#endif // INDIGOMONOMERADAPTER_H
