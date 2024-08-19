# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
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
