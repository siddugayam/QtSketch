#include "BondStereoPerception.h"
#include "indigo.h"
#include <QByteArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QString>
#include <algorithm>

// Indigo's JSON saver gives each disconnected component its own "mol0"/"mol1"/... node
// with LOCAL atom indices. The node-local numbers are translated back to `mol`'s own
// atom indices below, so the result is valid for multi-component molecules too.
std::map<std::pair<int,int>, QChar> computeIndigoBondCIP(int mol) {
    std::map<std::pair<int,int>, QChar> result;
    // Assumes this runs in a dedicated/throwaway Indigo session (matching
    // IndigoService.cpp's identical pattern) -- this option is process/session-global
    // and is never reset, so setting it on the shared main session would make every
    // future indigoJson()/toKetJson() call on that session also emit "cip" fields.
    indigoSetOptionBool("json-saving-add-stereo-desc", 1);
    const char* ketStr = indigoJson(mol);
    if (!ketStr) return result;
    QJsonDocument ketDoc = QJsonDocument::fromJson(QByteArray(ketStr));
    QJsonObject ketRoot = ketDoc.object();
    // The KET writer emits one "molN" node per connected component, and the atom numbers
    // inside each node are LOCAL to that component (0..k-1, in the component's own atom
    // order). The map returned here is keyed by the atom indices of `mol` itself, so each
    // local number has to be translated back through the component's atom list. Taking the
    // local numbers as they are was only right for a single-fragment molecule; with two or
    // more fragments (a salt, a reaction, two structures on one canvas) every label after
    // the first fragment landed on an unrelated atom pair or was lost.
    QHash<int, QList<int>> componentAtoms;   // component index -> its atoms' indices in `mol`, ascending
    {
        int iter = indigoIterateAtoms(mol);
        if (iter >= 0) {
            int a;
            while ((a = indigoNext(iter)) > 0) {
                componentAtoms[indigoComponentIndex(a)].append(indigoIndex(a));
                indigoFree(a);
            }
            indigoFree(iter);
        }
    }
    QJsonArray nodes = ketRoot.value("root").toObject().value("nodes").toArray();
    int moleculeOrdinal = 0;
    for (const QJsonValue &nodeVal : nodes) {
        QString ref = nodeVal.toObject().value("$ref").toString();
        if (ref.isEmpty()) continue;
        QJsonObject molObj = ketRoot.value(ref).toObject();
        if (molObj.value("type").toString() != "molecule") continue;
        const QList<int> atomsOfThisNode = componentAtoms.value(moleculeOrdinal++);
        // Guard against a writer that splits differently: only trust the translation when
        // the node holds exactly this component's atoms.
        if (molObj.value("atoms").toArray().size() != atomsOfThisNode.size()) continue;
        const QJsonArray bonds = molObj.value("bonds").toArray();
        for (const QJsonValue &bondVal : bonds) {
            QJsonObject bondObj = bondVal.toObject();
            QString label = bondObj.value("cip").toString();
            if (label != "E" && label != "Z") continue;
            QJsonArray bondAtoms = bondObj.value("atoms").toArray();
            if (bondAtoms.size() != 2) continue;
            const int l1 = bondAtoms.at(0).toInt();
            const int l2 = bondAtoms.at(1).toInt();
            if (l1 < 0 || l2 < 0 || l1 >= atomsOfThisNode.size() || l2 >= atomsOfThisNode.size()) continue;
            int a1 = atomsOfThisNode.at(l1);
            int a2 = atomsOfThisNode.at(l2);
            result[{std::min(a1, a2), std::max(a1, a2)}] = label == "Z" ? QChar('Z') : QChar('E');
        }
    }
    return result;
}
