#pragma once

#include <QString>

// SystemLocant: Structured interior and peripheral locant for general IUPAC / CAS rendering
// (FR-5.5.2 / P-25.3.3.1.2-.3).
// Used for generated names of ortho- and peri-fused systems.
struct SystemLocant {
    int peripheralNumber{0};
    char letter{'\0'};
    int superscript{0};
    int casPeripheralNumber{0};
    char casLetter{'\0'};

    bool operator==(const SystemLocant &other) const {
        return peripheralNumber == other.peripheralNumber &&
               letter == other.letter &&
               superscript == other.superscript;
    }

    bool operator!=(const SystemLocant &other) const {
        return !(*this == other);
    }

    bool operator<(const SystemLocant &other) const {
        if (peripheralNumber != other.peripheralNumber)
            return peripheralNumber < other.peripheralNumber;
        if (letter != other.letter)
            return letter < other.letter;
        return superscript < other.superscript;
    }
};

inline QString toSuperscriptDigits(int val) {
    if (val <= 0) return QString();
    QString s = QString::number(val);
    QString res;
    for (QChar c : s) {
        switch (c.toLatin1()) {
            case '0': res += QString::fromUtf8("\u2070"); break;
            case '1': res += QString::fromUtf8("\u00B9"); break;
            case '2': res += QString::fromUtf8("\u00B2"); break;
            case '3': res += QString::fromUtf8("\u00B3"); break;
            case '4': res += QString::fromUtf8("\u2074"); break;
            case '5': res += QString::fromUtf8("\u2075"); break;
            case '6': res += QString::fromUtf8("\u2076"); break;
            case '7': res += QString::fromUtf8("\u2077"); break;
            case '8': res += QString::fromUtf8("\u2078"); break;
            case '9': res += QString::fromUtf8("\u2079"); break;
            default: break;
        }
    }
    return res;
}

inline QString renderIupac(const SystemLocant &loc) {
    QString res = QString::number(loc.peripheralNumber);
    if (loc.letter != '\0') {
        res += QChar(loc.letter);
    }
    if (loc.superscript > 0) {
        res += toSuperscriptDigits(loc.superscript);
    }
    return res;
}

inline QString renderCas(const SystemLocant &loc) {
    if (loc.casPeripheralNumber != 0) {
        QString res = QString::number(loc.casPeripheralNumber);
        if (loc.casLetter != '\0') {
            res += QChar(loc.casLetter);
        }
        return res;
    }
    QString res = QString::number(loc.peripheralNumber);
    if (loc.letter != '\0') {
        res += QChar(loc.letter);
    }
    return res;
}
