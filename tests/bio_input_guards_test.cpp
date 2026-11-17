// Checks the IDT and HELM input guards that keep crash-prone text away from Indigo.
#include "BioInputGuards.h"
#include <cstdio>

static int g_pass = 0, g_fail = 0;
#define CHECK(cond) do { if (cond) ++g_pass; else { ++g_fail; std::printf("FAIL %s:%d %s\n", __FILE__, __LINE__, #cond); } } while (0)

static void test_idt() {
    QString bad;
    // Accepted: plain bases, inosine, RNA/2'-OMe/LNA prefixes, phosphorothioate, modifications.
    CHECK(BioInputGuards::idtIsSafe("ACGT", bad));
    CHECK(BioInputGuards::idtIsSafe("ACIT", bad));
    CHECK(BioInputGuards::idtIsSafe("rArCrG", bad));
    CHECK(BioInputGuards::idtIsSafe("mA+CA*C*G", bad));
    CHECK(BioInputGuards::idtIsSafe("/5Phos/ACGT/3BioTEG/", bad));
    // Rejected: ambiguity codes, which crash the macOS Indigo build.
    bad.clear();
    CHECK(!BioInputGuards::idtIsSafe("ACNT", bad));
    CHECK(bad == "N");
    bad.clear();
    CHECK(!BioInputGuards::idtIsSafe("/5Phos/ARYT", bad));
    CHECK(bad == "R,Y");
    CHECK(!BioInputGuards::idtIsSafe("mN", bad));
    CHECK(!BioInputGuards::idtIsSafe("AC(N)T", bad));
}

static void test_helm() {
    CHECK(BioInputGuards::helmIsSafe("PEPTIDE1{A.G.C}$$$$V2.0"));
    CHECK(BioInputGuards::helmIsSafe("RNA1{R(A)P.R(C)P}$$$$V2.0"));
    // Connection section commas and '+' inside inline SMILES are fine.
    CHECK(BioInputGuards::helmIsSafe("PEPTIDE1{A.C}|PEPTIDE2{C.A}$PEPTIDE1,PEPTIDE2,2:R3-1:R3$$$V2.0"));
    CHECK(BioInputGuards::helmIsSafe("PEPTIDE1{[C[N+](C)(C)C].A}$$$$V2.0"));
    CHECK(BioInputGuards::helmIsSafe("PEPTIDE1{A:0.5+G}$$$$V2.0"));
    // Parenthesised mixtures and alternatives crash Indigo.
    CHECK(!BioInputGuards::helmIsSafe("RNA1{R(A,C)P}$$$$V2.0"));
    CHECK(!BioInputGuards::helmIsSafe("RNA1{R(A+C)P}$$$$V2.0"));
    CHECK(!BioInputGuards::helmIsSafe("PEPTIDE1{(A,G).C}$$$$V2.0"));
}

int main() {
    test_idt();
    test_helm();
    std::printf("bio_input_guards_test: %d passed, %d failed\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
