#include "IupacNamer.h"
#include "NamingTrace.h"
#include "indigo.h"
#include "indigo-inchi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <QString>

static std::string escapeJson(const std::string &s) {
    std::string result;
    for (char c : s) {
        if (c == '"') result += "\\\"";
        else if (c == '\\') result += "\\\\";
        else if (c == '\n') result += "\\n";
        else if (c == '\t') result += "\\t";
        else result += c;
    }
    return result;
}

class FileTraceSink : public NamingTraceSink {
    FILE *m_fp;
public:
    explicit FileTraceSink(FILE *fp) : m_fp(fp) {}
    ~FileTraceSink() override { if (m_fp && m_fp != stdout && m_fp != stderr) fclose(m_fp); }
    void record(const NamingTraceEvent &event) override {
        std::string cat = event.category.toStdString();
        std::string det = event.detail.toStdString();
        fprintf(m_fp, "{\"category\":\"%s\",\"detail\":\"%s\"}\n",
            escapeJson(cat).c_str(), escapeJson(det).c_str());
        fflush(m_fp);
    }
};

int main(int argc, char** argv) {
    if (argc > 1 && strcmp(argv[1], "--batch") == 0) {
        if (argc < 3) {
            std::cerr << "usage: smiles_driver --batch <tsv_file>" << std::endl;
            return 2;
        }
        std::ifstream infile(argv[2]);
        if (!infile.is_open()) {
            std::cerr << "ERROR: cannot open batch file: " << argv[2] << std::endl;
            return 1;
        }
        qword sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        indigoInchiInit(sid);
        std::string line;
        int processedCount = 0;
        while (std::getline(infile, line)) {
            if (line.empty()) continue;
            processedCount++;
            if (processedCount % 200 == 0) {
                indigoReleaseSessionId(sid);
                sid = indigoAllocSessionId();
                indigoSetSessionId(sid);
                indigoInchiInit(sid);
            }
            std::stringstream ss(line);
            std::string chembl_id, drug_name, smiles;
            if (std::getline(ss, chembl_id, '\t') &&
                std::getline(ss, drug_name, '\t') &&
                std::getline(ss, smiles, '\t')) {
                int mol = indigoLoadMoleculeFromString(smiles.c_str());
                if (mol < 0) {
                    std::cout << chembl_id << "\t" << drug_name << "\t" << smiles << "\tLOAD_FAILED\t0\t" << indigoGetLastError() << "\n";
                    continue;
                }
                int heavyAtoms = 0;
                int iter = indigoIterateAtoms(mol);
                if (iter >= 0) {
                    int a;
                    while ((a = indigoNext(iter)) != 0) {
                        if (a == -1) break;
                        heavyAtoms++;
                        indigoFree(a);
                    }
                    indigoFree(iter);
                }
                IupacResult res = IupacNamer::generateName(mol);
                indigoFree(mol);
                if (res.success) {
                    std::cout << chembl_id << "\t" << drug_name << "\t" << smiles << "\tNAME\t" << heavyAtoms << "\t" << res.name.toStdString() << "\n";
                } else {
                    std::cout << chembl_id << "\t" << drug_name << "\t" << smiles << "\tREJECTED\t" << heavyAtoms << "\t" << res.error.toStdString() << "\n";
                }
            }
        }
        indigoReleaseSessionId(sid);
        return 0;
    }

    const char *tracePath = nullptr;
    int argOffset = 1;
    
    if (argc > 1 && strcmp(argv[1], "--trace-json") == 0) {
        if (argc < 3) {
            std::cerr << "usage: smiles_driver [--trace-json <path>] <smiles>" << std::endl;
            return 2;
        }
        tracePath = argv[2];
        argOffset = 3;
    }
    
    if (argc - argOffset < 1) {
        std::cerr << "usage: smiles_driver [--trace-json <path>] <smiles>" << std::endl;
        return 2;
    }
    
    FileTraceSink *sink = nullptr;
    
    if (tracePath) {
        FILE *traceFp = fopen(tracePath, "w");
        if (!traceFp) {
            std::cerr << "ERROR: cannot open trace file: " << tracePath << std::endl;
            return 1;
        }
        sink = new FileTraceSink(traceFp);
        setNamingTraceSink(sink);
    }
    
    qword sid = indigoAllocSessionId();
    indigoSetSessionId(sid);
    indigoInchiInit(sid);
    int mol = indigoLoadMoleculeFromString(argv[argOffset]);
    if (mol < 0) {
        std::cout << "LOAD_FAILED\t" << indigoGetLastError() << std::endl;
        if (sink) { delete sink; sink = nullptr; setNamingTraceSink(nullptr); }
        return 1;
    }
    int heavyAtoms = 0;
    int iter = indigoIterateAtoms(mol);
    if (iter >= 0) {
        int a;
        while ((a = indigoNext(iter)) != 0) {
            if (a == -1) break;
            heavyAtoms++;
            indigoFree(a);
        }
        indigoFree(iter);
    }
    IupacResult res = IupacNamer::generateName(mol);
    if (res.success) {
        std::cout << "NAME\t" << heavyAtoms << "\t" << res.name.toStdString() << std::endl;
    } else {
        std::cout << "REJECTED\t" << heavyAtoms << "\t" << res.error.toStdString() << std::endl;
    }
    
    if (sink) {
        delete sink;
        sink = nullptr;
        setNamingTraceSink(nullptr);
    }
    
    return 0;
}
