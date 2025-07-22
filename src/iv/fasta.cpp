// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <clice/clice.h>
#include <fstream>
#include <ivio/ivio.h>

namespace {
auto cli = clice::Argument {
    .args   = "fasta",
    .desc   = "process fasta files",
};


namespace merge {
void app_merge();
auto cliMerge = clice::Argument {
    .parent = &cli,
    .args   = "merge",
    .desc   = "merges fasta files",
    .value  = std::vector<std::filesystem::path>{},
    .cb     = app_merge,
};
auto cliOutput = clice::Argument {
    .parent = &cliMerge,
    .args   = {"-o", "--output"},
    .desc   = "path to a fasta file",
    .value  = std::string{"/dev/stdout"},
};
void app_merge() {
    auto writer = ivio::fasta::writer{{.output = *cliOutput}};
    size_t i{0};
    for (auto inputFile : *cliMerge) {
        ++i;
        fmt::print(stderr, "processed file {} of {} files\n", i, cliMerge->size());
        auto reader = ivio::fasta::reader{{.input = inputFile}};
        for (auto record : reader) {
            writer.write(record);
        }
    }
}
}


namespace filter {
void app_filter();
auto cliFilter = clice::Argument {
    .parent = &cli,
    .args   = "filter",
    .desc   = "filters or transforms fasta files",
    .value  = std::vector<std::filesystem::path>{},
    .cb     = app_filter,
};
auto cliOutput = clice::Argument {
    .parent = &cliFilter,
    .args   = {"-o", "--output"},
    .desc   = "path to a fasta file",
    .value  = std::string{"/dev/stdout"},
};
auto cliUpperCase = clice::Argument {
    .parent = &cliFilter,
    .args   = {"-u", "--upper-case"},
    .desc   = "translates words to upper case",
};
auto cliRandomize = clice::Argument {
    .parent = &cliFilter,
    .args   = {"-r", "--randomize"},
    .desc   = "a list of letters, if they not appear, replace with one of them randomly",
    .value  = std::string{},
};
auto cliTruncate = clice::Argument {
    .parent = &cliFilter,
    .args   = {"-t", "--truncate"},
    .desc   = "truncate after given numbers of nucleotides",
    .value  = std::numeric_limits<size_t>::max(),
};
void app_filter() {
    auto writer = ivio::fasta::writer{{.output = *cliOutput}};
    size_t i{0};
    auto bitmap = std::array<bool, 256>{};
    for (auto c : *cliRandomize) {
        bitmap[(uint8_t)c] = true;
    }
    if (cliRandomize && cliRandomize->size() == 0) {
        throw std::runtime_error{"-r, --random can not be an empty string"};
    }

    size_t count{};
    for (auto inputFile : *cliFilter) {
        ++i;
        fmt::print(stderr, "processed file {} of {} files\n", i, cliFilter->size());
        auto reader = ivio::fasta::reader{{.input = inputFile}};
        for (auto rec : reader) {
            auto r = ivio::fasta::record(rec);
            if (cliUpperCase) {
                for (auto& c : r.seq) {
                    c = std::toupper(c);
                }
            }
            if (cliRandomize) {
                for (auto& _c : r.seq) {
                    auto c = static_cast<uint8_t>(_c);
                    if (!bitmap[c]) {
                        auto i = rand() % cliRandomize->size();
                        _c = cliRandomize->at(i);
                    }
                }
            }
            if (count + r.seq.size() > *cliTruncate) {
                r.seq = r.seq.substr(0, *cliTruncate - count);
            }
            writer.write(r);
            count += r.seq.size();
            if (count == *cliTruncate) break;
        }
    }
}
}
}
