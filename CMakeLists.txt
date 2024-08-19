# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.16)

project (IVtools VERSION 1.0.0)

# Force language standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(OpenMP QUIET)
include(cmake/CPM.cmake)

find_package(fmt REQUIRED)

set(BUILD_DOC OFF)
set(BUILD_SANDBOX OFF)
set(SKIP_PERFORMANCE_COMPARISON ON)
CPMAddPackage("gh:USCiLab/cereal@1.3.2")
set(LIBSAIS_USE_OPENMP ON)
CPMAddPackage("gh:SGSSGene/libsais#8a5a468b2f31afe1db81888302ce4532c1436ece@2.8.4-1")

find_package(fmindex_collection REQUIRED PATHS lib/fmindex-collection)
find_package(clice REQUIRED PATHS lib/clice)
find_package(ivio REQUIRED PATHS lib/IVio)
find_package(ivsigma REQUIRED PATHS lib/IVSigma)

add_subdirectory(src/iv)
