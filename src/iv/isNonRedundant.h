// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <search_schemes/isComplete.h>

namespace search_schemes {

namespace complete::detail {

inline auto uniqueCover(Scheme const& ss, ErrorConfig const& errorConfig) -> size_t {
    size_t ct{};
    for (auto const& s : ss) {
        if (covers(s, errorConfig)) {
            ct += 1;
        }
    }
    return ct == 1;
}

/* checks if Scheme covers patterns up to a certain error
 *
 */
inline auto isUnique(Scheme const& ss, size_t minK, size_t maxK) -> bool {
    if (ss.empty()) return false;
    bool unique{true};
    auto len = ss.at(0).pi.size();
    generateErrorConfig([&](ErrorConfig const& config) {
        bool r = uniqueCover(ss, config);
        unique = unique && r;
    }, len, minK, maxK);
    return unique;
}
}

inline auto isNonRedundant(Scheme const& ss, size_t minK, size_t maxK) -> bool {
    return complete::detail::isUnique(ss, minK, maxK);
}

}
