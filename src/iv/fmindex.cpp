// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <clice/clice.h>
#include <fstream>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <libsais64.h>
#include <fmindex-collection/fmindex/FMIndex.h>
#include <fmindex-collection/fmindex/merge.h>
#include <fmt/std.h>

#include <unistd.h>

namespace {
void app();
auto cli = clice::Argument {
    .args   = "fmindex",
    .desc   = "construct and print an fmindex from fasta file",
    .value  = std::filesystem::path{},
    .cb     = app,
};

auto cliOutput = clice::Argument {
    .parent = &cli,
    .args   = {"-o", "--output"},
    .desc   = "path to a text file",
    .value  = std::filesystem::path{},
    .tags   = {"required"},
};

auto loadFastaFile(std::filesystem::path input) -> std::vector<std::vector<uint8_t>> {
    auto reader = ivio::fasta::reader{{.input = input}};
    auto res = std::vector<std::vector<uint8_t>>{};
    for (auto record_view : reader) {
        auto output = std::vector<uint8_t>{};
        for (auto c : record_view.seq) {
            output.push_back(c);
        }
        res.emplace_back(std::move(output));
    }
    return res;
}

void app() {
    auto seqs = loadFastaFile(*cli);

    fmt::print("start constructing\n");
    auto fullIndex = fmc::FMIndex</*.Sigma=*/256>{seqs, /*.samplingRate=*/1, /*.threads=*/1};
    fmt::print("writing to {}\n", *cliOutput);

    auto ofs = std::ofstream{*cliOutput};
    for (size_t i{0}; i < fullIndex.size(); ++i) {
        auto [sa, offset] = fullIndex.locate(i);
        auto [id, pos] = sa;
        auto _pos = pos;
        pos += offset;
        auto at = [&](size_t i) -> char {
            if (i == seqs[id].size()) return '$';
            auto j = i % (1+seqs[id].size());
            return (char)seqs[id][j];
        };
        auto substr = [&](size_t b, size_t e) -> std::string {
            auto r = std::string{};
            for (size_t i{b}; i < e; ++i) {
                r += at(i);
            }
            return r;
        };
        auto firstLetter = at(pos);
        auto midPart = substr(pos+1, pos+seqs[id].size());
//        auto midPart = str.substr(pos+1, pos+seqs[id].size()-2);
        auto lastLetter  = at(pos+seqs[id].size());
        ofs << fmt::format("{: >3}/{: >3}   {} {} {}\n", id, _pos, firstLetter, midPart, lastLetter);
    }
    ofs.close();
    fmt::print("finished\n");
}
}
