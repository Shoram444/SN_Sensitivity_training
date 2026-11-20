#include "read_generated_counts.C"

#include <cstdlib>

namespace {
std::string csv_field(const std::string& value) {
    std::string result = "\"";
    for (char c : value) {
        if (c == '"') result += '"';
        result += c;
    }
    return result + '"';
}
}

void list_root_events() {
    try {
        const char* list_path = std::getenv("ROOT_EVENTS_LIST");
        const char* output_path = std::getenv("ROOT_EVENTS_OUTPUT");
        if (!list_path || !output_path) throw std::runtime_error("missing input or output path");

        std::ifstream input(list_path, std::ios::binary);
        std::ofstream output(output_path);
        if (!input || !output) throw std::runtime_error("cannot open input list or output CSV");
        output << "file_path,nb_events\n";

        std::string path;
        Long64_t files = 0, failures = 0;
        while (std::getline(input, path, '\0')) {
            while (path.compare(0, 2, "./") == 0) path.erase(0, 2);
            ++files;
            output << csv_field(path) << ',';
            try {
                std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
                if (!file || file->IsZombie()) throw std::runtime_error("cannot open ROOT file");
                output << generated_count(file->Get("last_event_number"));
            } catch (const std::exception& error) {
                std::cerr << path << ": " << error.what() << '\n';
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
        std::cerr << "ROOT event listing failed: " << error.what() << '\n';
        gSystem->Exit(1);
    }
}
