// src/app/molecule/MoleculeSnapshot.h
#ifndef MOLECULESNAPSHOT_H
#define MOLECULESNAPSHOT_H

// One undo/redo entry: a full indigoClone of the molecule handle plus a
// value-copy of ExtensionData and the ID tables. Move-only; frees its clone
// on destruction. Must NOT outlive the EditableMolecule that created it
// (the clone lives on that molecule's Indigo session).

#include "ExtensionData.h"
#include <QHash>
#include <QMap>
#include <QString>

class EditableMolecule;

class MoleculeSnapshot {
public:
    MoleculeSnapshot(MoleculeSnapshot&& other) noexcept { *this = std::move(other); }
    MoleculeSnapshot& operator=(MoleculeSnapshot&& other) noexcept {
        if (this != &other) {
            freeHandle();
            m_clone = other.m_clone; m_session = other.m_session;
            m_ext = std::move(other.m_ext);
            m_atomIdx = std::move(other.m_atomIdx);
            m_bondIdx = std::move(other.m_bondIdx);
            m_nextAtomId = other.m_nextAtomId; m_nextBondId = other.m_nextBondId;
            m_sgroupIdx = std::move(other.m_sgroupIdx);
            m_sgroupExpanded = std::move(other.m_sgroupExpanded);
            m_nextSGroupId = other.m_nextSGroupId;
            m_sgroupLabels = std::move(other.m_sgroupLabels);
            m_ezWasSpecified = std::move(other.m_ezWasSpecified);
            other.m_clone = -1;
        }
        return *this;
    }
    MoleculeSnapshot(const MoleculeSnapshot&) = delete;
    MoleculeSnapshot& operator=(const MoleculeSnapshot&) = delete;
    ~MoleculeSnapshot() { freeHandle(); }

    bool isValid() const { return m_clone >= 0; }

private:
    friend class EditableMolecule;
    MoleculeSnapshot() = default;
    void freeHandle(); // indigoSetSessionId(m_session); indigoFree(m_clone) -- in .cpp side via EditableMolecule.cpp include

    int m_clone = -1;
    unsigned long long m_session = 0;
    ExtensionData m_ext;
    QHash<int, int> m_atomIdx;
    QHash<int, int> m_bondIdx;
    int m_nextAtomId = 1;
    int m_nextBondId = 1;
    QHash<int, int> m_sgroupIdx;
    QHash<int, bool> m_sgroupExpanded;
    int m_nextSGroupId = 1;
    QHash<int, QString> m_sgroupLabels;
    QMap<int, int> m_ezWasSpecified;
};

#endif // MOLECULESNAPSHOT_H
