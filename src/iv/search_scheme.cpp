// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include <clice/clice.h>
#include <search_schemes/search_schemes.h>

#include "error_fmt.h"

namespace {
void app();
auto cli                 = clice::Argument{ .args   = "search_scheme",
                                            .desc   = "generates and info about search schemes",
                                            .cb     = app,
};
auto cliListGenerator    = clice::Argument{ .parent = &cli,
                                            .args   = {"list-generators"},
                                            .desc   = "show a list of generators"
};
auto cliGenerator        = clice::Argument{ .parent = &cli,
                                            .args   = {"-g", "--generator"},
                                            .desc   = "which generator to use?",
                                            .value  = std::string{"pigeon"},
};
auto cliQueryLength      = clice::Argument{ .parent = &cli,
                                            .args   = {"-l", "--length"},
                                            .desc   = "the assumed query length, when applying node count",
                                            .value  = int64_t{150},
};
auto cliReferenceLength  = clice::Argument{ .parent = &cli,
                                            .args   = {"--ref-length"},
                                            .desc   = "the assumed length of the refernce text",
                                            .value  = int64_t{1'000'000'000},
};
auto cliMinAllowedErrors = clice::Argument{ .parent = &cli,
                                            .args   = {"--min-error"},
                                            .desc   = "minimum errors that have to appear, such that the search scheme accepts it",
                                            .value  = int64_t{0},
};
auto cliMaxAllowedErrors = clice::Argument{ .parent = &cli,
                                            .args   = {"-k", "--max-error"},
                                            .desc   = "maximum errors that can appear",
                                            .value  = int64_t{2},
};
auto cliAlphabetSize     = clice::Argument{ .parent = &cli,
                                            .args   = {"--sigma"},
                                            .desc   = "Size of the alphabet, e.g.: '4' for ACGT or  '5' for ACGTN",
                                            .value  = int64_t{4},
};
auto cliAll              = clice::Argument{ .parent = &cli,
                                            .args   = {"-a", "--all"},
                                            .desc   = "print information table about all generators",
};
auto cliYaml             = clice::Argument{ .parent = &cli,
                                            .args   = {"-y", "--yaml"},
                                            .desc   = "print in a yaml compatible format",
};


void printSingleScheme() {
    // pick the correct generator
    auto iter = search_schemes::generator::all.find(*cliGenerator);
    if (iter == search_schemes::generator::all.end()) {
        throw error_fmt{"can not find generator \"{}\"", *cliGenerator};
    }
    auto const& e = iter->second;

    // generate search schemes
    auto sss = e.generator(*cliMinAllowedErrors, *cliMaxAllowedErrors, *cliAlphabetSize, *cliReferenceLength);

    // expand search schemes, so they match the length of the query
    auto ss = expand(sss, *cliQueryLength);

    // expand search schemes, but "dynamically" trying to minimize weighted node count
    auto dss = expandByWNC</*Edit=*/true>(sss, *cliQueryLength, *cliAlphabetSize, *cliReferenceLength);

    auto parts = sss[0].pi.size();

    fmt::print("# Search Scheme Information\n");
    fmt::print("name:                {}\n", e.name);
    fmt::print("description:         {}\n", e.description);
    fmt::print("alphabet size:       {}\n", *cliAlphabetSize);
    fmt::print("min errors:          {}\n", *cliMinAllowedErrors);
    fmt::print("max errors:          {}\n", *cliMaxAllowedErrors);
    fmt::print("reference length:    {}\n", *cliReferenceLength);
    fmt::print("number of parts:     {}\n", parts);
    fmt::print("number of searches:  {}\n", ss.size());
    fmt::print("valid:               {}\n", isValid(sss));
    fmt::print("complete:            {}\n", isComplete(sss, *cliMinAllowedErrors, *cliMaxAllowedErrors));
    fmt::print("node count:          {}\n", nodeCount</*Edit=*/false>(ss, *cliAlphabetSize));
    fmt::print("weighted node count: {}\n", weightedNodeCount</*Edit=*/false>(ss, *cliAlphabetSize, *cliReferenceLength));
    fmt::print("dynamic wnc:         {}\n", weightedNodeCount</*Edit=*/false>(dss, *cliAlphabetSize, *cliReferenceLength));

    fmt::print("searches:  {:^{}}  {:^{}}  {:^{}}\n", "pi", parts*3, "L", parts*3, "U", parts*3);
    for (auto const& s : sss) {
        fmt::print("           {{{}}}, {{{}}}, {{{}}}\n", fmt::join(s.pi, ", "), fmt::join(s.l, ", "), fmt::join(s.u, ", "));
    }
}

void printTable() {
    fmt::print("# Search Scheme Information\n");
    fmt::print("alphabet size:       {}\n", *cliAlphabetSize);
    fmt::print("min errors:          {}\n", *cliMinAllowedErrors);
    fmt::print("max errors:          {}\n", *cliMaxAllowedErrors);
    fmt::print("reference length:    {}\n", *cliReferenceLength);

    fmt::print("{:^15} | {:^6} {:^8} {:^6} {:^8} | {:^30} | {:^24}\n", "name", "parts", "searches", "valid", "complete", "node count ham/edit", "weighted node count ham/edit");
    auto order = std::vector<std::string>{"backtracking", "optimum", "01*0", "01*0_opt", "pigeon", "pigeon_opt", "suffix", "h2-k1", "h2-k2", "h2-k3", "kianfar", "kucherov-k1", "kucherov-k2"};
    for (auto const& [key, e] : search_schemes::generator::all) {
        if (std::find(order.begin(), order.end(), key) == order.end()) {
            order.push_back(key);
            fmt::print("WARNING: missing {} in order list\n", key);
        }
    }
    for (auto const& o : order) {
        if (auto iter = search_schemes::generator::all.find(o); iter == search_schemes::generator::all.end()) {
            fmt::print("Warning: generator {} doesn't exists\n", o);
            continue;
        }

        auto sigma = *cliAlphabetSize;
        auto N     = *cliReferenceLength;

        auto const& e = search_schemes::generator::all[o];
        // generate search schemes
        auto sss = e.generator(*cliMinAllowedErrors, *cliMaxAllowedErrors, sigma, N);

        // expand search schemes, so they match the length of the query
        auto ss = expand(sss, *cliQueryLength);

        // expand by node count
        auto dss_ham  = expandByNC</*Edit=*/false>(sss, *cliQueryLength, sigma);
        auto dss_edit = expandByNC</*Edit=*/true> (sss, *cliQueryLength, sigma);

        // expand by weighted node count
        auto dess_ham  = expandByWNC</*Edit=*/false>(sss, *cliQueryLength, sigma, N);
        auto dess_edit = expandByWNC</*Edit=*/true> (sss, *cliQueryLength, sigma, N);

        auto parts = ss.size()>0?sss[0].pi.size():0;

        struct {
            long double countHam;
            long double countEdit;
        } stat_ss, stat_ss_w, stat_dss, stat_dess;

        auto complete = isComplete(sss, *cliMinAllowedErrors, *cliMaxAllowedErrors);
        auto valid    = isValid(sss);
        stat_ss   = { .countHam  = nodeCount</*Edit=*/false>(ss, sigma),
                      .countEdit = nodeCount</*Edit=*/true>(ss, sigma)};
        stat_ss_w = { .countHam  = weightedNodeCount</*Edit=*/false>(ss, sigma, N),
                      .countEdit = weightedNodeCount</*Edit=*/true>(ss, sigma, N)};
        stat_dss  = { .countHam  = nodeCount</*Edit=*/false>(dss_ham, sigma),
                      .countEdit = nodeCount</*Edit=*/true>(dss_edit, sigma)};
        stat_dess = { .countHam  = weightedNodeCount</*Edit=*/false>(dss_ham, sigma, N),
                      .countEdit = weightedNodeCount</*Edit=*/true>(dss_edit, sigma, N)};

        fmt::print("{:>15} | {:>6} {:>8} {:^6} {:^8} | {:>15.0f} {:>15.0f}  | {:>12.2f} {:>12.2f} | {:>15.0f} {:>15.0f} | {:>12.2f} {:>12.2f}\n", e.name, parts, sss.size(), valid, complete, stat_ss.countHam, stat_ss.countEdit, stat_ss_w.countHam, stat_ss_w.countEdit, stat_dss.countHam, stat_dss.countEdit, stat_dess.countHam, stat_dess.countEdit);
    }
}
void printYaml() {
    fmt::print("# Search Scheme Information\n");
    fmt::print("alphabet size:       {}\n", *cliAlphabetSize);
    fmt::print("min errors:          {}\n", *cliMinAllowedErrors);
    fmt::print("max errors:          {}\n", *cliMaxAllowedErrors);
    fmt::print("reference length:    {}\n", *cliReferenceLength);
    fmt::print("---\n");
    for (auto k{*cliMinAllowedErrors}; k <= *cliMaxAllowedErrors; ++k) {
        for (auto const& [key, e] : search_schemes::generator::all) {
            // generate search schemes
            auto sss = e.generator(*cliMinAllowedErrors, k, *cliAlphabetSize, *cliReferenceLength);

            // expand search schemes, so they match the length of the query
            auto ss = expand(sss, *cliQueryLength);

            auto name          = e.name;
            auto parts         = ss.size()>0?sss[0].pi.size():0;
            auto searches      = ss.size();
            auto valid         = isValid(sss);
            auto complete      = isComplete(sss, *cliMinAllowedErrors, k);
            auto nodeCount     = search_schemes::nodeCount</*Edit=*/false>(ss, *cliAlphabetSize);
            auto weightedCount = weightedNodeCount</*Edit=*/false>(ss, *cliAlphabetSize, *cliReferenceLength);
            fmt::print("- name: \"{}\"\n", name);
            fmt::print("  parts: {}\n", parts);
            fmt::print("  searchCt: {}\n", searches);
            fmt::print("  valid: {}\n", valid);
            fmt::print("  complete: {}\n", complete);
            fmt::print("  nodeCount: {}\n", nodeCount);
            fmt::print("  weightedNodeCount: {:.2f}\n", weightedCount);
            fmt::print("  searches:\n");
            for (auto const& s : sss) {
                fmt::print("  - pi: [{}]\n", fmt::join(s.pi, ", "));
                fmt::print("    l: [{}]\n", fmt::join(s.l, ", "));
                fmt::print("    u: [{}]\n", fmt::join(s.u, ", "));
            }
        }
    }
}




void app() {
    if (cliListGenerator) {
        for (auto const& [key, e] : search_schemes::generator::all) {
            fmt::print("{:>15} - {}\n", e.name, e.description);
        }
        return;
    }

    if (cliAll && cliYaml) {
        printYaml();
    } else if (cliAll) {
        printTable();
    } else {
        printSingleScheme();
    }

}
}
