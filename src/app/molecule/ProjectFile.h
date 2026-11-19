// src/app/molecule/ProjectFile.h
#ifndef PROJECTFILE_H
#define PROJECTFILE_H

#include <QString>
#include <QJsonObject>
#include <QJsonArray>
#include <QJsonDocument>
#include "DrawingStyle.h"
#include "PageSettings.h"

// Versioned .qtsketch project file format
// 
// Format specification:
// Version 1 (legacy):
// {
//   "formatVersion": 1,
//   "chemistry": {"molfile": "..."},
//   "drawingStyle": { ... DrawingStyle fields ... },
//   "pageSettings": { ... PageSettings fields ... },
//   "modelToMmScale": 3.386666...,
//   "ezIntent": { bondId1: stereoValue, ... },
//   "sgroups": [ {id, expanded, label, atoms}, ... ],
//   "annotations": [ {type, id, x, y, content, bold, italic}, ... ],
//   "images": [ {id, x, y, w, h, dataBase64}, ... ],
//   "rxnArrows": [ {id, x1, y1, x2, y2, mode, above, below}, ... ],
//   "rxnPluses": [ {id, x, y}, ... ],
//   "multitailArrows": [ {id, points: [x1,y1, ...]}, ... ],
//   "rGroups": [ {number, logic: {range, resth, ifthen}, fragments}, ... ],
//   "brackets": [ {minX, minY, maxX, maxY}, ... ]
// }
//
// Version 2 (with biopolymer support):
// Adds:
//   "biopolymer": {  // Optional, only present if document has biopolymer data
//     "chains": [
//       {
//         "id": 1,
//         "name": "A",
//         "polymerType": "PEPTIDE",
//         "nTerminus": { ... },
//         "cTerminus": { ... },
//         "residues": [
//           {
//             "id": 1,
//             "symbol": "A",
//             "alias": "Ala",
//             "sequencePosition": 1,
//             "atomIds": [1, 2, 3, ...]
//           },
//           ...
//         ]
//       },
//       ...
//     ],
//     "connections": [
//       {
//         "id": 1,
//         "fromResidueId": 1,
//         "toResidueId": 2,
//         "bondId": 1,
//         "type": "peptide"
//       },
//       ...
//     ],
//     "nextIds": {
//       "residue": 10,
//       "chain": 2,
//       "connection": 5
//     }
//   }

class EditableMolecule;
class BiopolymerState;

class ProjectFile {
public:
    // Current format version
    static constexpr int CurrentFormatVersion = 2;
    
    // Minimum supported version for loading
    static constexpr int MinimumSupportedVersion = 1;
    
    // Format version 1 (legacy)
    static constexpr int LegacyFormatVersion = 1;

    // Save a document to .qtsketch format
    static bool saveProject(const QString& filePath, 
                            const EditableMolecule& molecule,
                            const DrawingStyle& drawingStyle,
                            const PageSettings& pageSettings);

    // Load a document from .qtsketch format
    // Returns empty string on success, error message on failure
    static QString loadProject(const QString& filePath,
                               EditableMolecule& molecule,
                               DrawingStyle& drawingStyle,
                               PageSettings& pageSettings);
    
    // Load a document from .qtsketch format with biopolymer state
    // Returns empty string on success, error message on failure
    static QString loadProjectWithBiopolymer(const QString& filePath,
                                           EditableMolecule& molecule,
                                           DrawingStyle& drawingStyle,
                                           PageSettings& pageSettings,
                                           BiopolymerState& biopolymerState);

    // Serialize document to JSON (for autosave/recovery)
    static QJsonObject toJson(const EditableMolecule& molecule,
                              const DrawingStyle& drawingStyle,
                              const PageSettings& pageSettings);
    
    // Serialize document to JSON with biopolymer state
    static QJsonObject toJsonWithBiopolymer(const EditableMolecule& molecule,
                                            const DrawingStyle& drawingStyle,
                                            const PageSettings& pageSettings,
                                            const BiopolymerState& biopolymerState);

    // Deserialize document from JSON
    // Returns true on success, false on failure
    static bool fromJson(const QJsonObject& json,
                         EditableMolecule& molecule,
                         DrawingStyle& drawingStyle,
                         PageSettings& pageSettings);
    
    // Deserialize document from JSON with biopolymer state
    // Returns true on success, false on failure
    static bool fromJsonWithBiopolymer(const QJsonObject& json,
                                       EditableMolecule& molecule,
                                       DrawingStyle& drawingStyle,
                                       PageSettings& pageSettings,
                                       BiopolymerState& biopolymerState);

    // Check if a file is a valid .qtsketch file
    static bool isValidProjectFile(const QString& filePath);

    // Check if a file is a legacy .ket autosave file
    static bool isLegacyRecoveryFile(const QString& filePath);

    // Get format version from a file
    static int getFormatVersion(const QString& filePath);

    // Validate that a JSON object contains a valid project
    static bool validateProjectJson(const QJsonObject& json, QString& errorMessage);
};

#endif // PROJECTFILE_H
