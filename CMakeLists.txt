# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.16)

project (IVtools LANGUAGES CXX C
                 VERSION 1.0.0)

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(OpenMP QUIET)
include(cmake/CPM.cmake)
CPMAddPackage("gh:SGSSGene/cpm.dependencies@1.0.0")
CPMLoadDependenciesFile("${CMAKE_CURRENT_SOURCE_DIR}/cpm.dependencies")

add_subdirectory(src/iv)
