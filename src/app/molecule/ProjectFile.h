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

class EditableMolecule;

class ProjectFile {
public:
    // Current format version
    static constexpr int CurrentFormatVersion = 1;
    
    // Minimum supported version for loading
    static constexpr int MinimumSupportedVersion = 1;

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

    // Serialize document to JSON (for autosave/recovery)
    static QJsonObject toJson(const EditableMolecule& molecule,
                              const DrawingStyle& drawingStyle,
                              const PageSettings& pageSettings);

    // Deserialize document from JSON
    // Returns true on success, false on failure
    static bool fromJson(const QJsonObject& json,
                         EditableMolecule& molecule,
                         DrawingStyle& drawingStyle,
                         PageSettings& pageSettings);

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
