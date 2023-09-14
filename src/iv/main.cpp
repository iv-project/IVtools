// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include <clice/clice.h>
#include <fmt/format.h>

namespace {
auto cliHelp = clice::Argument { .args     = {"-h", "--help"},
                                 .desc     = "prints the help page",
                                 .cb       = []{ fmt::print("{}", clice::generateHelp()); exit(0); },
};
}

void slix_script_main(std::string script);

int main(int argc, char** argv) {
    try {
        if (auto failed = clice::parse(argc, argv); failed) {
            fmt::print(stderr, "parsing failed {}\n", *failed);
            return 1;
        }
    } catch (std::exception const& e) {
        fmt::print(stderr, "error {}\n", e.what());
    }
    return 0;
}
