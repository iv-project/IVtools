// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <clice/clice.h>
#include <fstream>
#include <ivio/ivio.h>

namespace {
auto cli = clice::Argument{ .args   = "fasta",
                            .desc   = "merge fasta files",
};


void app_merge();
auto cliMerge = clice::Argument{ .parent = &cli,
                                 .args   = "merge",
                                 .value  = std::vector<std::filesystem::path>{},
                                 .cb     = app_merge,
};
auto cliOutput = clice::Argument{ .parent = &cliMerge,
                                  .args   = {"-o", "--output"},
                                  .desc   = "path to a fasta file",
                                  .value  = std::string{"/dev/stdout"},
};
void app_merge() {
    auto writer = ivio::fasta::writer{{.output = *cliOutput}};
    size_t i{0};
    for (auto inputFile : *cliMerge) {
        ++i;
        fmt::print("processed file {} of {} files\n", i, cliMerge->size());
        auto reader = ivio::fasta::reader{{.input = inputFile}};
        for (auto record : reader) {
            writer.write(record);
        }
    }
}
}
