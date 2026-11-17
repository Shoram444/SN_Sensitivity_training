#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include "TFile.h"
#include "TParameter.h"
#include "TSystem.h"

namespace {
Long64_t generated_count(TObject* object) {
    if (!object) throw std::runtime_error("missing last_event_number");
    Long64_t count;
    if (auto* p = dynamic_cast<TParameter<Long64_t>*>(object)) {
        count = p->GetVal();
    } else if (auto* p = dynamic_cast<TParameter<Int_t>*>(object)) {
        count = p->GetVal();
    } else if (auto* p = dynamic_cast<TParameter<Long_t>*>(object)) {
        count = p->GetVal();
    } else if (auto* p = dynamic_cast<TParameter<Double_t>*>(object)) {
        const double value = p->GetVal();
        if (!std::isfinite(value) || value <= 0 || value != std::floor(value) ||
            value > 9007199254740991.0) {
            throw std::runtime_error("last_event_number is not an exact positive integer");
        }
        count = static_cast<Long64_t>(value);
    } else {
        throw std::runtime_error(std::string("unsupported last_event_number type: ") + object->ClassName());
    }
    if (count <= 0) throw std::runtime_error("last_event_number must be positive");
    return count;
}
}

void read_generated_counts(const char* list_path, const char* output_path) {
    try {
        std::ifstream input(list_path);
        std::ofstream output(output_path);
        if (!input || !output) throw std::runtime_error("cannot open metadata input/output list");
        output << "file_index\tlast_event_number\n";
        std::string path;
        Long64_t index = 0;
        while (std::getline(input, path)) {
            if (path.empty()) throw std::runtime_error("empty ROOT path");
            std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
            if (!file || file->IsZombie()) throw std::runtime_error("cannot open ROOT file: " + path);
            try {
                output << index << '\t' << generated_count(file->Get("last_event_number")) << '\n';
            } catch (const std::exception& error) {
                throw std::runtime_error(path + ": " + error.what());
            }
            if (!output) throw std::runtime_error("failed writing metadata report");
            ++index;
        }
        if (input.bad()) throw std::runtime_error("failed reading ROOT path list");
        output.close();
        if (!output) throw std::runtime_error("failed closing metadata report");
        std::cout << "Read generated totals from " << index << " files; no event trees read.\n";
    } catch (const std::exception& error) {
        std::cerr << "Metadata extraction failed: " << error.what() << '\n';
        gSystem->Exit(1);
    }
}
