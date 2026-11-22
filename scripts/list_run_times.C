#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>
#include "TFile.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TTree.h"

namespace {
std::string csv_field(const std::string& value) {
    std::string result = "\"";
    for (char c : value) {
        if (c == '"') result += '"';
        result += c;
    }
    return result + '"';
}

double run_time(TObject* object) {
    if (!object) throw std::runtime_error("missing run_time");
    double value;
    if (auto* p = dynamic_cast<TParameter<Double_t>*>(object)) {
        value = p->GetVal();
    } else if (auto* p = dynamic_cast<TParameter<Float_t>*>(object)) {
        value = p->GetVal();
    } else if (auto* p = dynamic_cast<TParameter<Long64_t>*>(object)) {
        value = p->GetVal();
    } else if (auto* p = dynamic_cast<TParameter<Long_t>*>(object)) {
        value = p->GetVal();
    } else if (auto* p = dynamic_cast<TParameter<Int_t>*>(object)) {
        value = p->GetVal();
    } else {
        throw std::runtime_error(std::string("unsupported run_time type: ") + object->ClassName());
    }
    if (!std::isfinite(value) || value < 0) throw std::runtime_error("run_time is not a finite non-negative number");
    return value;
}

// The first event's run_number; empty when the tree has no events.
std::string tree_run_number(TFile& file) {
    auto* tree = dynamic_cast<TTree*>(file.Get("Event"));
    if (!tree || tree->GetEntries() == 0) return "";
    TLeaf* leaf = tree->GetLeaf("run_number");
    if (!leaf) return "";
    leaf->GetBranch()->GetEntry(0);
    return std::to_string(static_cast<Long64_t>(leaf->GetValue()));
}

// Run number from the parent directory, .../run_XXX/extracted_data.root.
std::string directory_run_number(const std::string& path) {
    static const std::regex pattern("run_0*([0-9]+)");
    const std::string parent = path.substr(0, path.find_last_of('/'));
    const std::string name = parent.substr(parent.find_last_of('/') + 1);
    std::smatch match;
    if (!std::regex_match(name, match, pattern)) throw std::runtime_error("parent directory is not run_XXX");
    return match[1].str();
}
}

void list_run_times() {
    try {
        const char* list_path = std::getenv("RUN_TIMES_LIST");
        const char* output_path = std::getenv("RUN_TIMES_OUTPUT");
        if (!list_path || !output_path) throw std::runtime_error("missing input or output path");

        std::ifstream input(list_path, std::ios::binary);
        std::ofstream output(output_path);
        if (!input || !output) throw std::runtime_error("cannot open input list or output CSV");
        output << "file_path,run_number,run_time\n" << std::setprecision(17);

        std::string path;
        Long64_t files = 0, failures = 0;
        while (std::getline(input, path, '\0')) {
            ++files;
            output << csv_field(path) << ',';
            try {
                std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
                if (!file || file->IsZombie()) throw std::runtime_error("cannot open ROOT file");
                const std::string run = directory_run_number(path);
                const double time = run_time(file->Get("run_time"));
                const std::string from_tree = tree_run_number(*file);
                if (!from_tree.empty() && from_tree != run) {
                    std::cerr << path << ": Event run_number " << from_tree << " differs from directory\n";
                }
                output << run << ',' << time;
            } catch (const std::exception& error) {
                std::cerr << path << ": " << error.what() << '\n';
                output << ',';
                ++failures;
            }
            output << '\n';
            if (!output) throw std::runtime_error("failed writing CSV");
        }
        if (input.bad()) throw std::runtime_error("failed reading ROOT path list");
        if (files == 0) throw std::runtime_error("no ROOT files found");
        output.close();
        if (!output) throw std::runtime_error("failed closing CSV");
        std::cerr << "Read " << files << " files; " << failures << " failed.\n";
        if (failures != 0) gSystem->Exit(1);
    } catch (const std::exception& error) {
        std::cerr << "Run time listing failed: " << error.what() << '\n';
        gSystem->Exit(1);
    }
}
