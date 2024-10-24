#include <clice/clice.h>
#include <fstream>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <libsais64.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/fmindex/FMIndex.h>
#include <fmindex-collection/fmindex/merge.h>

#include <unistd.h>
auto getTotalSystemMemory() -> size_t {
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

namespace {
void app();
auto cli = clice::Argument{ .args   = "bwt",
                            .desc   = "construct a bwt from a fasta file",
                            .cb     = app,
};

auto cliInput = clice::Argument{ .parent = &cli,
                                 .args   = {"-i", "--input"},
                                 .desc   = "path to a fasta file",
                                 .value  = std::filesystem::path{},
                                 .tags   = {"required"},
};

auto cliOutput = clice::Argument{ .parent = &cli,
                                  .args   = {"-o", "--output"},
                                  .desc   = "path to a text file",
                                  .value  = std::string{},
                                  .tags   = {"required"},
};


auto cliMaxMemoryUsage = clice::Argument{ .parent = &cli,
                                          .args   = {"-m", "--max_memory"},
                                          .desc   = "maximum memory usage for suffix array",
                                          .value  = getTotalSystemMemory()/2,
};

auto cliAsciiMode      = clice::Argument{ .parent = &cli,
                                          .args   = {"--ascii"},
                                          .desc   = "assumes a text file",
};


char randomPick() {
    switch(rand()% 4) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
    }
    throw std::runtime_error("should never happen");
}

auto loadFastaFile(std::filesystem::path input) -> std::vector<std::vector<uint8_t>> {
    auto reader = ivio::fasta::reader{{.input = input}};
    auto res = std::vector<std::vector<uint8_t>>{};
    for (auto record_view : reader) {
        auto output = ivs::convert_char_to_rank<ivs::d_dna5>(record_view.seq);
        for (auto& c : output) {
            if (!ivs::verify_rank(c)) {
                c = 5; // 'N'
            }
        }
        res.emplace_back(std::move(output));
    }
    return res;
}

void app_dna5() {
    fmt::print("reading string T from fasta file...\n");
    auto seqs = loadFastaFile(*cliInput);

    #if FMC_USE_SDSL
    using OccTable = fmindex_collection::occtable::Sdsl_wt_bldc</*.Sigma=*/256>;
    #else
    using OccTable = fmindex_collection::occtable::EprV7</*.Sigma=*/256>;
    #endif
    using Index    = fmindex_collection::FMIndex<OccTable>;

    auto fullIndex = Index{};

    fmt::print("start constructing\n");
    size_t processedSeqs{};
    while (processedSeqs < seqs.size()) {
        auto next = std::vector<std::vector<uint8_t>>{};
        size_t acc{};
        auto requiredMem = (acc + seqs[processedSeqs].size())*8;
        while (next.empty() || requiredMem < *cliMaxMemoryUsage) {
            if (requiredMem > *cliMaxMemoryUsage) {
                fmt::print("Using more memory than allowed {} -> {}\n", *cliMaxMemoryUsage, requiredMem);
            }
            next.push_back(std::move(seqs[processedSeqs]));
            acc += next.back().size();
            processedSeqs += 1;
            if (processedSeqs == seqs.size()) break;
            requiredMem = (acc + seqs[processedSeqs].size())*8;
        }

        fmt::print("processing {} of {} sequences ({})\n", processedSeqs, seqs.size(), acc+next.size());
        auto partIndex = Index{next, /*.samplingRate =*/ 8, /*.threadNbr =*/ 16};

        fmt::print("merging... ({} + {})\n", fullIndex.size(), partIndex.size());
        if (fullIndex.size() > 0) {
            fullIndex = fmindex_collection::merge(fullIndex, partIndex);
        } else {
            fullIndex = std::move(partIndex);
        }
    }
    fmt::print("writing to {}\n", *cliOutput);
    auto ofs = std::ofstream{*cliOutput};
    for (size_t i{0}; i < fullIndex.size(); ++i) {
        auto r = fullIndex.occ.symbol(i);
        auto c = ivs::d_dna5::rank_to_char(r);
        ofs << c;
    }
    ofs.close();
    fmt::print("finished\n");
}

auto loadTextFile(std::filesystem::path input) -> std::vector<std::vector<uint8_t>> {
    auto res = std::vector<std::vector<uint8_t>>{};
    auto line = std::string{};
    auto fs = std::ifstream{input};
    while (std::getline(fs, line)) {
        auto o = std::vector<uint8_t>{};
        o.resize(line.size());
        std::memcpy(o.data(), line.data(), o.size());
        res.push_back(std::move(o));
    }
    return res;
}


void app_ascii() {
    fmt::print("reading string T from text file...\n");
    auto seqs = loadTextFile(*cliInput);

    #if FMC_USE_SDSL
    using OccTable = fmindex_collection::occtable::Sdsl_wt_bldc</*.Sigma=*/256>;
    #else
    using OccTable = fmindex_collection::occtable::EprV7</*.Sigma=*/256>;
    #endif
    using Index    = fmindex_collection::FMIndex<OccTable>;

    auto fullIndex = Index{};

    fmt::print("start constructing\n");
    size_t processedSeqs{};
    while (processedSeqs < seqs.size()) {
        auto next = std::vector<std::vector<uint8_t>>{};
        size_t acc{};
        auto requiredMem = (acc + seqs[processedSeqs].size())*8;
        while (next.empty() || requiredMem < *cliMaxMemoryUsage) {
            if (requiredMem > *cliMaxMemoryUsage) {
                fmt::print("Using more memory than allowed {} -> {}\n", *cliMaxMemoryUsage, requiredMem);
            }
            next.push_back(std::move(seqs[processedSeqs]));
            acc += next.back().size();
            processedSeqs += 1;
            if (processedSeqs == seqs.size()) break;
            requiredMem = (acc + seqs[processedSeqs].size())*8;
        }

        fmt::print("processing {} of {} sequences ({})\n", processedSeqs, seqs.size(), acc+next.size());
        auto partIndex = Index{next, /*.samplingRate =*/ 8, /*.threadNbr =*/ 16};

        fmt::print("merging... ({} + {})\n", fullIndex.size(), partIndex.size());
        if (fullIndex.size() > 0) {
            fullIndex = fmindex_collection::merge(fullIndex, partIndex);
        } else {
            fullIndex = std::move(partIndex);
        }
    }
    fmt::print("writing to {}\n", *cliOutput);
    auto ofs = std::ofstream{*cliOutput};
    for (size_t i{0}; i < fullIndex.size(); ++i) {
        auto r = fullIndex.occ.symbol(i);
        ofs.write((char*)&r, 1);
    }
    ofs.close();
    fmt::print("finished\n");
}


void app() {
    if (cliAsciiMode) {
        app_ascii();
    } else {
        app_dna5();
    }
}

}
