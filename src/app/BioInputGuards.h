#pragma once
// Input checks run before IDT and HELM text reaches Indigo. The macOS Indigo
// build segfaults on:
//   IDT:  IUPAC ambiguity bases R Y K S W H B V D N, e.g. "ACNT" (M loads
//         empty instead, so it is rejected too)
//   HELM: a parenthesised mixture or alternative such as "(A,G)" or "R(A+C)P"
// Both are rejected here with a readable message instead.

#include <QSet>
#include <QString>
#include <QStringList>

namespace BioInputGuards {

// Returns false and lists the offending codes in badChars if the IDT text holds
// an ambiguity base outside a /modification/ block.
inline bool idtIsSafe(const QString &text, QString &badChars) {
    static const QSet<QChar> ambiguity = {
        'R','Y','K','M','S','W','H','B','V','D','N'
    };
    QSet<QChar> bad;
    bool inMod = false;
    for (QChar ch : text) {
        if (ch == QLatin1Char('/')) { inMod = !inMod; continue; }
        if (!inMod && ambiguity.contains(ch)) bad.insert(ch);
    }
    if (bad.isEmpty()) return true;
    QStringList list;
    for (QChar c : bad) list << QString(c);
    list.sort();
    badChars = list.join(QLatin1Char(','));
    return false;
}

// Returns false if a HELM simple-polymer block {...} holds a parenthesised
// group containing ',' or '+'. Inline SMILES [..] and quoted text are skipped.
inline bool helmIsSafe(const QString &text) {
    int braceDepth = 0, parenDepth = 0, bracketDepth = 0;
    bool inQuote = false;
    for (QChar ch : text) {
        if (inQuote) { if (ch == QLatin1Char('"')) inQuote = false; continue; }
        if (ch == QLatin1Char('"')) { inQuote = true; continue; }
        if (ch == QLatin1Char('[')) { ++bracketDepth; continue; }
        if (ch == QLatin1Char(']')) { if (bracketDepth > 0) --bracketDepth; continue; }
        if (bracketDepth > 0) continue;
        if (ch == QLatin1Char('{')) { ++braceDepth; continue; }
        if (ch == QLatin1Char('}')) { if (braceDepth > 0) --braceDepth; parenDepth = 0; continue; }
        if (braceDepth == 0) continue;
        if (ch == QLatin1Char('(')) { ++parenDepth; continue; }
        if (ch == QLatin1Char(')')) { if (parenDepth > 0) --parenDepth; continue; }
        if (parenDepth > 0 && (ch == QLatin1Char(',') || ch == QLatin1Char('+'))) return false;
    }
    return true;
}

} // namespace BioInputGuards
