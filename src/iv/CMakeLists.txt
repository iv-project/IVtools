# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

cmake_minimum_required (VERSION 3.16)

add_executable(iv
    columba_prepare.cpp
    inspect_fasta.cpp
    main.cpp
    read_simulator.cpp
    search_scheme.cpp
    bwt.cpp
    fasta.cpp
)

target_link_libraries(iv
    fmindex-collection::fmindex-collection
    fmt::fmt
    ivio::ivio
    ivsigma::ivsigma
    clice::clice
)

target_compile_features(iv PUBLIC cxx_std_20)
