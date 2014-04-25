// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/detail/has_identity.hpp>
#include <boost/units_blas/operations.hpp>

#include <boost/test/minimal.hpp>


struct derived_from_E_matrix_2x2_double_type :
    E_matrix_2x2_double_type
{};

struct derived_from_E_matrix_2x2_time_type :
    E_matrix_2x2_time_type
{};

typedef bub::matrix<
    std::tuple<time_, time_>,
    std::tuple<time_, double>
> uninvertible_matrix_2x2;

int test_main (int, char *[])
{
    static_assert(bub::detail::has_identity<E_matrix_2x2_double_type>::value, "");
    static_assert(bub::detail::has_identity<derived_from_E_matrix_2x2_double_type>::value, "");
    static_assert(bub::detail::has_identity<E_matrix_2x2_time_type>::value, "");
    static_assert(bub::detail::has_identity<derived_from_E_matrix_2x2_time_type>::value, "");

    static_assert(!bub::detail::has_identity<int>::value, "");
    static_assert(!bub::detail::has_identity<uninvertible_matrix_2x2>::value, "");

    return 0;
}
