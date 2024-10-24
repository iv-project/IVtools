#include <clice/clice.h>
#include <fstream>
#include <ivio/ivio.h>

namespace {
void app_merge();
auto cli = clice::Argument{ .args   = "fasta",
                            .desc   = "merge fasta files",
};

auto cliMerge = clice::Argument{ .parent = &cli,
                                 .args   = "merge",
                                 .cb     = app_merge,
};

auto cliInput = clice::Argument{ .parent = &cliMerge,
                                 .args   = {"-i", "--input"},
                                 .desc   = "path to multiple fasta files",
                                 .value  = std::vector<std::filesystem::path>{},
                                 .tags   = {"required"},
};

auto cliOutput = clice::Argument{ .parent = &cliMerge,
                                  .args   = {"-o", "--output"},
                                  .desc   = "path to a fasta file",
                                  .value  = std::string{},
                                  .tags   = {"required"},
};

void app_merge() {
    auto writer = ivio::fasta::writer{{.output = *cliOutput}};
    size_t i{0};
    for (auto inputFile : *cliInput) {
        ++i;
        fmt::print("processed file {} of {} files\n", i, cliInput->size());
        auto reader = ivio::fasta::reader{{.input = inputFile}};
        for (auto record : reader) {
            writer.write(record);
        }
    }
}
}
