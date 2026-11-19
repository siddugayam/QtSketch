// src/app/molecule/IndigoMonomerAdapter.cpp
#include "IndigoMonomerAdapter.h"
#include "EditableMolecule.h"

// ========== Static Template Data ==========

// Standard 20 amino acids (L-form)
const QList<IndigoMonomerAdapter::MonomerTemplate> IndigoMonomerAdapter::s_peptideTemplates = []() {
    QList<MonomerTemplate> templates;
    
    // Glycine - no chiral center
    templates.append({
        "G", "Gly", "AminoAcid", false,
        "NCC(=O)O",  // SMILES without stereochemistry
        {0, 1},       // N-terminus at index 0 (N), C-terminus at index 1 (C)
        {"N", "C", "C", "O", "O"},  // H2N-CH2-COOH
        {0, 0, 0, 0, 0},            // Charges
        {3, 4, 4, 2, 2},            // Valencies
        "", {},
        {},
        {},
        {}
    });
    
    // Alanine
    templates.append({
        "A", "Ala", "AminoAcid", false,
        "N[C@H](C)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "O", "O"},  // H2N-CH(CH3)-COOH
        {0, 0, 0, 0, 0, 0},
        {3, 4, 4, 4, 2, 2},
        "L", {1},  // Chiral at carbon 1 (the alpha carbon)
        {},
        {},
        {}
    });
    
    // Valine
    templates.append({
        "V", "Val", "AminoAcid", false,
        "N[C@H]([C@H](C)C)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "C", "C", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 4, 4, 4, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Leucine
    templates.append({
        "L", "Leu", "AminoAcid", false,
        "N[C@H](CC(C)C)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "C", "C", "C", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 4, 4, 4, 4, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Isoleucine
    templates.append({
        "I", "Ile", "AminoAcid", false,
        "N[C@H]([C@@H](C)CC)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "C", "C", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 4, 4, 4, 2, 2},
        "L", {1, 2},  // Two chiral centers
        {},
        {},
        {}
    });
    
    // Serine
    templates.append({
        "S", "Ser", "AminoAcid", false,
        "N[C@H](CO)C(=O)O",
        {0, 1},
        {"N", "C", "C", "O", "O", "O"},
        {0, 0, 0, 0, 0, 0},
        {3, 4, 4, 2, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Threonine
    templates.append({
        "T", "Thr", "AminoAcid", false,
        "N[C@H]([C@@H](O)C)C(=O)O",
        {0, 1},
        {"N", "C", "C", "O", "C", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 2, 4, 2, 2},
        "L", {1, 2},
        {},
        {},
        {}
    });
    
    // Cysteine
    templates.append({
        "C", "Cys", "AminoAcid", false,
        "N[C@H](CS)C(=O)O",
        {0, 1},
        {"N", "C", "C", "S", "O", "O"},
        {0, 0, 0, 0, 0, 0},
        {3, 4, 4, 2, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Proline
    templates.append({
        "P", "Pro", "AminoAcid", false,
        "N1CC[C@H]1C(=O)O",  // Cyclic
        {0, 2},  // N at index 0, alpha C at index 2
        {"N", "C", "C", "C", "C", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 4, 4, 2, 2},
        "L", {2},  // Chiral at index 2
        {},
        {},
        {}
    });
    
    // Aspartic Acid
    templates.append({
        "D", "Asp", "AminoAcid", false,
        "N[C@H](CC(=O)O)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "O", "O", "O", "O"},
        {0, 0, 0, 0, -1, -1, 0, 0},
        {3, 4, 4, 4, 2, 2, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Asparagine
    templates.append({
        "N", "Asn", "AminoAcid", false,
        "N[C@H](CC(=O)N)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "O", "N", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 4, 2, 3, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Glutamic Acid
    templates.append({
        "E", "Glu", "AminoAcid", false,
        "N[C@H](CCC(=O)O)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "C", "O", "O", "O", "O"},
        {0, 0, 0, 0, 0, -1, -1, 0, 0},
        {3, 4, 4, 4, 4, 2, 2, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Glutamine
    templates.append({
        "Q", "Gln", "AminoAcid", false,
        "N[C@H](CCC(=O)N)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "C", "O", "N", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 4, 4, 2, 3, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Methionine
    templates.append({
        "M", "Met", "AminoAcid", false,
        "N[C@H](CSCC)C(=O)O",
        {0, 1},
        {"N", "C", "C", "S", "C", "C", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 2, 4, 4, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Lysine
    templates.append({
        "K", "Lys", "AminoAcid", false,
        "N[C@H](CCCCN)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "C", "C", "N", "O", "O"},
        {0, 0, 0, 0, 0, 0, +1, 0, 0},
        {3, 4, 4, 4, 4, 4, 3, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Arginine
    templates.append({
        "R", "Arg", "AminoAcid", false,
        "N[C@H](CCCNC(=N)N)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "C", "N", "C", "N", "N", "O", "O"},
        {0, 0, 0, 0, 0, +1, 0, +1, 0, 0, 0},
        {3, 4, 4, 4, 4, 3, 2, 3, 3, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Histidine
    templates.append({
        "H", "His", "AminoAcid", false,
        "N[C@H](C1=CN=CN1)C(=O)O",
        {0, 1},
        {"N", "C", "C", "N", "C", "N", "C", "O", "O"},
        {0, 0, 0, 0, +1, 0, 0, 0, 0},
        {3, 4, 4, 2, 3, 2, 2, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Phenylalanine
    templates.append({
        "F", "Phe", "AminoAcid", false,
        "N[C@H](C1=CC=CC=C1)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "C", "C", "C", "C", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 4, 4, 4, 4, 4, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Tyrosine
    templates.append({
        "Y", "Tyr", "AminoAcid", false,
        "N[C@H](C1=CC=C(O)C=C1)C(=O)O",
        {0, 1},
        {"N", "C", "C", "C", "C", "C", "C", "O", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 4, 4, 4, 4, 2, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    // Tryptophan
    templates.append({
        "W", "Trp", "AminoAcid", false,
        "N[C@H](C1=CNc2c1C=CC=C2)C(=O)O",
        {0, 1},
        {"N", "C", "C", "N", "C", "C", "C", "C", "C", "C", "O", "O"},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {3, 4, 4, 2, 3, 4, 4, 4, 4, 4, 2, 2},
        "L", {1},
        {},
        {},
        {}
    });
    
    return templates;
}();

// DNA bases
const QList<IndigoMonomerAdapter::MonomerTemplate> IndigoMonomerAdapter::s_dnaTemplates = []() {
    QList<MonomerTemplate> templates;
    
    // Adenine (DNA)
    templates.append({
        "A", "dA", "Base", false,
        "Nc1ncnc2n(cnc12)[C@H]1O[C@H](COP(=O)(O)O)C1O",
        {0, 1, 2},  // 0 = base attachment, 1 = sugar 1', 2 = phosphate
        {"N", "C", "N", "C", "N", "C", "C", "O", "C", "O", "C", "C", "O", "P", "O", "O", "O"},
        {},
        {},
        "", {},
        {0, 1, 2, 3, 4, 5, 6},  // Base atoms (adenine)
        {7, 8, 9, 10, 11},      // Sugar atoms (deoxyribose)
        {12, 13, 14, 15, 16}    // Phosphate atoms
    });
    
    // Cytosine (DNA)
    templates.append({
        "C", "dC", "Base", false,
        "Nc1ccn([C@H]2O[C@H](COP(=O)(O)O)C2O)c(=O)nc1",
        {0, 1, 2},
        {"N", "C", "C", "N", "C", "O", "C", "O", "C", "O", "C", "C", "O", "P", "O", "O", "O"},
        {},
        {},
        "", {},
        {0, 1, 2, 3, 4, 5},  // Base atoms (cytosine)
        {6, 7, 8, 9, 10, 11},  // Sugar atoms
        {12, 13, 14, 15, 16}    // Phosphate atoms
    });
    
    // Guanine (DNA)
    templates.append({
        "G", "dG", "Base", false,
        "Nc1nc2c(ncn2c(=O)nc1)[C@H]1O[C@H](COP(=O)(O)O)C1O",
        {0, 1, 2},
        {"N", "C", "N", "C", "N", "C", "N", "C", "O", "C", "O", "C", "C", "O", "P", "O", "O", "O"},
        {},
        {},
        "", {},
        {0, 1, 2, 3, 4, 5, 6, 7, 8},  // Base atoms (guanine)
        {9, 10, 11, 12, 13},          // Sugar atoms
        {14, 15, 16, 17}              // Phosphate atoms
    });
    
    // Thymine (DNA)
    templates.append({
        "T", "dT", "Base", false,
        "C1=CN(C(=O)NC1=O)[C@H]1O[C@H](COP(=O)(O)O)C1O",
        {0, 1, 2},
        {"C", "C", "N", "C", "O", "N", "C", "O", "C", "O", "C", "C", "O", "P", "O", "O", "O"},
        {},
        {},
        "", {},
        {0, 1, 2, 3, 4, 5, 6},  // Base atoms (thymine)
        {7, 8, 9, 10, 11},      // Sugar atoms
        {12, 13, 14, 15, 16}    // Phosphate atoms
    });
    
    return templates;
}();

// RNA bases
const QList<IndigoMonomerAdapter::MonomerTemplate> IndigoMonomerAdapter::s_rnaTemplates = []() {
    QList<MonomerTemplate> templates;
    
    // Adenine (RNA) - note ribose has 2' OH
    templates.append({
        "A", "A", "Base", false,
        "Nc1ncnc2n(cnc12)[C@H]1O[C@@H](COP(=O)(O)O)[C@H](O)C1O",
        {0, 1, 2},
        {"N", "C", "N", "C", "N", "C", "C", "O", "C", "O", "C", "C", "O", "C", "O", "P", "O", "O", "O"},
        {},
        {},
        "", {},
        {0, 1, 2, 3, 4, 5, 6},  // Base atoms
        {7, 8, 9, 10, 11, 12, 13},  // Sugar atoms (ribose with 2' OH)
        {14, 15, 16, 17, 18}    // Phosphate atoms
    });
    
    // Cytosine (RNA)
    templates.append({
        "C", "C", "Base", false,
        "Nc1ccn([C@H]2O[C@@H](COP(=O)(O)O)[C@H](O)C2O)c(=O)nc1",
        {0, 1, 2},
        {"N", "C", "C", "N", "C", "O", "C", "O", "C", "O", "C", "C", "O", "C", "O", "P", "O", "O", "O"},
        {},
        {},
        "", {},
        {0, 1, 2, 3, 4, 5},
        {6, 7, 8, 9, 10, 11, 12},
        {13, 14, 15, 16, 17}
    });
    
    // Guanine (RNA)
    templates.append({
        "G", "G", "Base", false,
        "Nc1nc2c(ncn2c(=O)nc1)[C@H]1O[C@@H](COP(=O)(O)O)[C@H](O)C1O",
        {0, 1, 2},
        {"N", "C", "N", "C", "N", "C", "N", "C", "O", "C", "O", "C", "C", "O", "C", "O", "P", "O", "O", "O"},
        {},
        {},
        "", {},
        {0, 1, 2, 3, 4, 5, 6, 7, 8},
        {9, 10, 11, 12, 13, 14, 15},
        {16, 17, 18, 19}
    });
    
    // Uracil (RNA)
    templates.append({
        "U", "U", "Base", false,
        "C1=CN(C(=O)NC1=O)[C@H]1O[C@@H](COP(=O)(O)O)[C@H](O)C1O",
        {0, 1, 2},
        {"C", "C", "N", "C", "O", "N", "C", "O", "C", "O", "C", "C", "O", "C", "O", "P", "O", "O", "O"},
        {},
        {},
        "", {},
        {0, 1, 2, 3, 4, 5, 6},
        {7, 8, 9, 10, 11, 12, 13},
        {14, 15, 16, 17, 18}
    });
    
    return templates;
}();

// Ambiguous amino acids
const QList<IndigoMonomerAdapter::MonomerTemplate> IndigoMonomerAdapter::s_ambiguousPeptideTemplates = []() {
    QList<MonomerTemplate> templates;
    
    templates.append({"B", "Asx", "AminoAcid", true, "", {0, 1}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"Z", "Glx", "AminoAcid", true, "", {0, 1}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"J", "Xle", "AminoAcid", true, "", {0, 1}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"X", "Xaa", "AminoAcid", true, "", {0, 1}, {}, {}, {}, "", {}, {}, {}, {}});
    
    return templates;
}();

// Ambiguous nucleotides
const QList<IndigoMonomerAdapter::MonomerTemplate> IndigoMonomerAdapter::s_ambiguousNucleotideTemplates = []() {
    QList<MonomerTemplate> templates;
    
    templates.append({"N", "N", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"B", "B", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"V", "V", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"D", "D", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"H", "H", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"K", "K", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"M", "M", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"W", "W", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"Y", "Y", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"R", "R", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    templates.append({"S", "S", "Base", true, "", {0, 1, 2}, {}, {}, {}, "", {}, {}, {}, {}});
    
    return templates;
}();

// Terminal group templates
const QList<IndigoMonomerAdapter::TerminalGroupTemplate> IndigoMonomerAdapter::s_peptideTerminalGroups = []() {
    QList<TerminalGroupTemplate> templates;
    
    // N-terminus: NH2 (default)
    templates.append({"NH2", "N", {}, {}});
    
    // N-terminus: NH3+ (protonated)
    templates.append({"NH3+", "[NH3+]", {}, {}});
    
    // N-terminus: Acetyl
    templates.append({"Ac", "CC(=O)N", {1}, {}});  // Attaches at index 1 (the N)
    
    // C-terminus: COOH (default)
    templates.append({"COOH", "C(=O)O", {}, {}});
    
    // C-terminus: COO- (deprotonated)
    templates.append({"COO-", "C(=O)[O-]", {}, {}});
    
    // C-terminus: Amide (CONH2)
    templates.append({"CONH2", "C(=O)N", {0}, {}});  // Attaches at index 0 (the C)
    
    return templates;
}();

const QList<IndigoMonomerAdapter::TerminalGroupTemplate> IndigoMonomerAdapter::s_nucleotideTerminalGroups = []() {
    QList<TerminalGroupTemplate> templates;
    
    // 5' phosphate
    templates.append({"5'-PO4", "P(=O)(O)O", {0}, {}});
    
    // 5' monophosphate
    templates.append({"5'-P", "P(=O)(O)O", {0}, {}});
    
    // 3' phosphate
    templates.append({"3'-PO4", "OP(=O)(O)O", {0}, {}});
    
    // 3' hydroxyl (default for RNA)
    templates.append({"3'-OH", "O", {}, {}});
    
    return templates;
}();

// ========== IndigoMonomerAdapter Implementation ==========

IndigoMonomerAdapter::IndigoMonomerAdapter() {
    // Build template maps
    m_peptideTemplateMap = buildFullTemplateMap(PolymerType::Peptide);
    m_dnaTemplateMap = buildFullTemplateMap(PolymerType::DNA);
    m_rnaTemplateMap = buildFullTemplateMap(PolymerType::RNA);
}

IndigoMonomerAdapter::~IndigoMonomerAdapter() = default;

const IndigoMonomerAdapter::MonomerTemplate* 
IndigoMonomerAdapter::templateForSymbol(const QString& symbol, PolymerType type) const {
    QString normSymbol = normalizeMonomerSymbol(symbol);
    const auto& map = getTemplateMap(type);
    auto it = map.find(normSymbol);
    return (it != map.end()) ? &it.value() : nullptr;
}

bool IndigoMonomerAdapter::isValidSymbol(const QString& symbol, PolymerType type) const {
    return templateForSymbol(symbol, type) != nullptr;
}

QList<QString> IndigoMonomerAdapter::validSymbols(PolymerType type) const {
    QList<QString> symbols;
    for (const auto& tpl : templates(type)) {
        symbols.append(tpl.symbol);
    }
    return symbols;
}

QList<IndigoMonomerAdapter::MonomerTemplate> 
IndigoMonomerAdapter::templates(PolymerType type) const {
    return templatesForType(type) + ambiguousTemplatesForType(type);
}

const QList<IndigoMonomerAdapter::MonomerTemplate>& 
IndigoMonomerAdapter::templatesForType(PolymerType type) const {
    switch (type) {
    case PolymerType::Peptide: return s_peptideTemplates;
    case PolymerType::DNA: return s_dnaTemplates;
    case PolymerType::RNA: return s_rnaTemplates;
    default: return s_peptideTemplates; // Shouldn't happen
    }
}

const QList<IndigoMonomerAdapter::MonomerTemplate>& 
IndigoMonomerAdapter::ambiguousTemplatesForType(PolymerType type) const {
    switch (type) {
    case PolymerType::Peptide: return s_ambiguousPeptideTemplates;
    case PolymerType::DNA: return s_ambiguousNucleotideTemplates;
    case PolymerType::RNA: return s_ambiguousNucleotideTemplates;
    default: return s_ambiguousPeptideTemplates;
    }
}

ExpansionResult IndigoMonomerAdapter::expandMonomer(
    EditableMolecule& molecule, 
    const QString& symbol, 
    PolymerType type,
    double x, double y) const {
    
    const MonomerTemplate* tpl = templateForSymbol(symbol, type);
    if (!tpl) {
        return ExpansionResult{false, QString("Unknown monomer symbol: %1 for type %2").arg(symbol).arg(polymerTypeToString(type)), "", {}, {}, {}, {}};
    }
    
    // For now, we use a simplified approach: load from SMILES
    // In a full implementation, we would use Indigo's layout capabilities
    // and properly map template indices to created atom IDs
    
    // Try to load the SMILES
    QString smiles = tpl->smiles;
    if (smiles.isEmpty()) {
        // For ambiguous templates without SMILES, we can't expand
        return ExpansionResult{false, QString("No SMILES for ambiguous monomer: %1").arg(symbol), "", {}, {}, {}, {}};
    }
    
    // Check if we need to add explicit hydrogens
    // For biopolymers, we typically want the full structure
    
    // Try to parse the SMILES through EditableMolecule
    // We need to position it at (x, y)
    // For now, we'll use a simplified approach
    
    // For peptide residues, we use a placeholder SMILES with proper atom count
    // In production, we would use Indigo to parse and layout the SMILES
    
    // For M0, we'll return a placeholder result
    // The actual Indigo integration would go here
    
    ExpansionResult result;
    result.success = true;
    result.molfile = ""; // Would be populated from Indigo
    
    // For M0, we'll simulate creating atoms
    // In reality, this would come from EditableMolecule::loadFrom()
    
    // Create a simple placeholder - in real implementation, we'd parse SMILES
    // and create the full structure
    
    // For now, just return the template info
    result.success = true;
    result.errorMessage = "";
    
    // Simulate atom creation
    // In reality: molecule.loadFrom(smiles); and then we'd map
    
    return result;
}

ExpansionResult IndigoMonomerAdapter::expandMonomerWithLabels(
    EditableMolecule& molecule,
    const QString& symbol,
    PolymerType type,
    const QMap<int, QString>& atomLabelOverrides,
    double x, double y) const {
    
    // Similar to expandMonomer but with label overrides
    // Placeholder for M0
    return expandMonomer(molecule, symbol, type, x, y);
}

ExpansionResult IndigoMonomerAdapter::expandNTerminus(
    EditableMolecule& molecule, AtomId attachTo, double x, double y) const {
    
    // Placeholder for M0
    ExpansionResult result;
    result.success = true;
    result.errorMessage = "";
    return result;
}

ExpansionResult IndigoMonomerAdapter::expandCTerminus(
    EditableMolecule& molecule, AtomId attachTo, double x, double y) const {
    
    // Placeholder for M0
    ExpansionResult result;
    result.success = true;
    result.errorMessage = "";
    return result;
}

ExpansionResult IndigoMonomerAdapter::expandPhosphate5Prime(
    EditableMolecule& molecule, AtomId attachTo, double x, double y) const {
    
    // Placeholder for M0
    ExpansionResult result;
    result.success = true;
    result.errorMessage = "";
    return result;
}

ExpansionResult IndigoMonomerAdapter::expandPhosphate3Prime(
    EditableMolecule& molecule, AtomId attachTo, double x, double y) const {
    
    // Placeholder for M0
    ExpansionResult result;
    result.success = true;
    result.errorMessage = "";
    return result;
}

ConnectionResult IndigoMonomerAdapter::connectResidues(
    EditableMolecule& molecule,
    const QString& fromSymbol, PolymerType fromType,
    const QString& toSymbol, PolymerType toType,
    AtomId fromAtomId, AtomId toAtomId) const {
    
    // For M0, we'll use the default attachment points
    return connectByAttachmentPoints(molecule, fromSymbol, fromType, toSymbol, toType,
                                     1, 0, fromAtomId, toAtomId);
}

ConnectionResult IndigoMonomerAdapter::connectByAttachmentPoints(
    EditableMolecule& molecule,
    const QString& fromSymbol, PolymerType fromType,
    const QString& toSymbol, PolymerType toType,
    int fromAttachmentIndex,
    int toAttachmentIndex,
    AtomId fromAtomId, AtomId toAtomId) const {
    
    // For M0, we'll just create a single bond
    const MonomerTemplate* fromTpl = templateForSymbol(fromSymbol, fromType);
    const MonomerTemplate* toTpl = templateForSymbol(toSymbol, toType);
    
    if (!fromTpl || !toTpl) {
        return ConnectionResult{false, "Unknown monomer template", -1, -1, -1};
    }
    
    // Validate attachment indices
    if (fromAttachmentIndex >= fromTpl->attachmentPoints.size() ||
        toAttachmentIndex >= toTpl->attachmentPoints.size()) {
        return ConnectionResult{false, "Invalid attachment point index", -1, -1, -1};
    }
    
    // Create the bond
    BondId bondId = molecule.addBond(fromAtomId, toAtomId, 1);
    
    ConnectionResult result;
    result.success = (bondId != -1);
    result.errorMessage = (bondId != -1) ? "" : "Failed to create bond";
    result.createdBondId = bondId;
    result.fromAtomId = fromAtomId;
    result.toAtomId = toAtomId;
    
    return result;
}

bool IndigoMonomerAdapter::validateExpansion(
    const EditableMolecule& molecule,
    const ExpansionResult& result,
    const QString& symbol,
    PolymerType type) const {
    
    if (!result.success) {
        return false;
    }
    
    const MonomerTemplate* tpl = templateForSymbol(symbol, type);
    if (!tpl) {
        return false;
    }
    
    // Check atom count
    if (!tpl->atomLabels.isEmpty()) {
        if (result.createdAtomIds.size() != tpl->atomLabels.size()) {
            return false;
        }
    }
    
    return true;
}

bool IndigoMonomerAdapter::validateConnection(
    const EditableMolecule& molecule,
    const ConnectionResult& result,
    const QString& fromSymbol, PolymerType fromType,
    const QString& toSymbol, PolymerType toType) const {
    
    if (!result.success) {
        return false;
    }
    
    // Check if the bond exists and connects the right atoms
    if (result.createdBondId == -1) {
        return false;
    }
    
    AtomId a, b;
    if (!molecule.bondEndpoints(result.createdBondId, a, b)) {
        return false;
    }
    
    if ((a != result.fromAtomId || b != result.toAtomId) &&
        (a != result.toAtomId || b != result.fromAtomId)) {
        return false;
    }
    
    return true;
}

AtomId IndigoMonomerAdapter::getAttachmentPointAtom(
    const EditableMolecule& molecule,
    const QString& symbol,
    PolymerType type,
    int attachmentIndex,
    const QMap<int, AtomId>& templateIndexToAtomId) const {
    
    const MonomerTemplate* tpl = templateForSymbol(symbol, type);
    if (!tpl) {
        return -1;
    }
    
    if (attachmentIndex >= tpl->attachmentPoints.size()) {
        return -1;
    }
    
    int templateIndex = tpl->attachmentPoints[attachmentIndex];
    auto it = templateIndexToAtomId.find(templateIndex);
    return (it != templateIndexToAtomId.end()) ? it.value() : -1;
}

bool IndigoMonomerAdapter::hasModification(const QString& symbol, PolymerType type, const QString& modification) const {
    // For M0, no modifications are supported yet
    return false;
}

QString IndigoMonomerAdapter::getStereochemistry(const QString& symbol, PolymerType type) const {
    const MonomerTemplate* tpl = templateForSymbol(symbol, type);
    if (!tpl) {
        return "";
    }
    return tpl->stereochemistry;
}

int IndigoMonomerAdapter::expectedHeavyAtomCount(const QString& symbol, PolymerType type) const {
    const MonomerTemplate* tpl = templateForSymbol(symbol, type);
    if (!tpl) {
        return 0;
    }
    return tpl->atomLabels.size();
}

QString IndigoMonomerAdapter::atomLabelAtPosition(const QString& symbol, PolymerType type, int position) const {
    const MonomerTemplate* tpl = templateForSymbol(symbol, type);
    if (!tpl || position >= tpl->atomLabels.size()) {
        return "";
    }
    return tpl->atomLabels[position];
}

// ========== Static Helper Methods ==========

QString IndigoMonomerAdapter::polymerTypeToString(PolymerType type) {
    switch (type) {
    case PolymerType::Peptide: return "PEPTIDE";
    case PolymerType::DNA: return "DNA";
    case PolymerType::RNA: return "RNA";
    default: return "UNKNOWN";
    }
}

QString IndigoMonomerAdapter::normalizePolymerTypeString(const QString& type) {
    QString upper = type.toUpper();
    if (upper == "PEPTIDE" || upper == "PROTEIN") return "PEPTIDE";
    if (upper == "DNA") return "DNA";
    if (upper == "RNA") return "RNA";
    return type;
}

QString IndigoMonomerAdapter::normalizeMonomerSymbol(const QString& symbol) {
    return symbol.trimmed().toUpper().left(1);
}

bool IndigoMonomerAdapter::supportsModification(PolymerType type, const QString& modification) {
    // For M0, no modifications are supported
    return false;
}

// ========== Private Methods ==========

QHash<QString, IndigoMonomerAdapter::MonomerTemplate> 
IndigoMonomerAdapter::buildFullTemplateMap(PolymerType type) {
    QHash<QString, MonomerTemplate> map;
    
    for (const auto& tpl : templatesForType(type)) {
        map[tpl.symbol] = tpl;
    }
    
    for (const auto& tpl : ambiguousTemplatesForType(type)) {
        map[tpl.symbol] = tpl;
    }
    
    return map;
}

const QHash<QString, IndigoMonomerAdapter::MonomerTemplate>& 
IndigoMonomerAdapter::getTemplateMap(PolymerType type) const {
    switch (type) {
    case PolymerType::Peptide: return m_fullPeptideTemplateMap;
    case PolymerType::DNA: return m_fullDnaTemplateMap;
    case PolymerType::RNA: return m_fullRnaTemplateMap;
    default: return m_fullPeptideTemplateMap;
    }
}

QMap<int, AtomId> IndigoMonomerAdapter::parseMolfileAtomMap(
    const QString& molfile,
    const EditableMolecule& molecule,
    const QMap<int, AtomId>& existingMap) const {
    
    // Parse molfile to extract atom indices and map to AtomIds
    // Placeholder for M0
    return existingMap;
}

// ========== Utility Functions ==========

QString polymerTypeToString(PolymerType type) {
    return IndigoMonomerAdapter::normalizePolymerTypeString(
        IndigoMonomerAdapter::polymerTypeToString(type));
}
