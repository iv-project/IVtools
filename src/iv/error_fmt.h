// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <fmt/format.h>
#include <fmt/std.h>
#include <stdexcept>

struct error_fmt : std::exception {
    std::string message;

    template <typename... Args>
    error_fmt(fmt::format_string<Args...> s, Args&&... args) {
        message = fmt::format(s, std::forward<Args>(args)...);
    }

    auto what() const noexcept -> char const* override {
        return message.c_str();
    }
};
