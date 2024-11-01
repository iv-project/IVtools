// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <clice/clice.h>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <random>

namespace {
void app();
auto cli = clice::Argument{ .args   = "inspect_fasta",
                            .desc   = "simulates reads of a certain length",
                            .cb     = app,
};

auto cliInput = clice::Argument{ .parent = &cli,
                                 .args   = {"-i", "--input"},
                                 .desc   = "path to a fasta file",
                                 .value  = std::filesystem::path{},
                                 .tags   = {"required"},
};

auto cliSeqId = clice::Argument{ .parent = &cli,
                                 .args   = {"-s", "--seq_id"},
                                 .desc   = "id of the sequence (as a number)",
                                 .value  = size_t{},
};
auto cliSeqPos = clice::Argument{ .parent = &cli,
                                  .args   = {"-p", "--seq_pos"},
                                  .desc   = "position inside the sequence",
                                  .value  = size_t{},
};

auto cliReadLength = clice::Argument{ .parent = &cli,
                                      .args   = {"-l", "--read_length"},
                                      .desc   = "length of the sequence to show",
                                      .value  = size_t{150},
};

void app() {
    size_t seqId{};
    for (auto record : ivio::fasta::reader {{*cliInput}}) {
        ++seqId;
        if (seqId < *cliSeqId+1) continue;
        auto seq = record.seq.substr(*cliSeqPos, *cliReadLength);
        fmt::print(">{}\n{}\n", record.id, seq);

        break;
    }
}
}
