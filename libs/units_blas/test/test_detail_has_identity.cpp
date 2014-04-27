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

using uninvertible_matrix_2x2 = bub::matrix<
    std::tuple<time_, time_>,
    std::tuple<time_, double>
>;

using uninvertible_matrix_3x3 = bub::matrix<
    std::tuple<length_sq,          length_sq,          length_sq_per_time>,
    std::tuple<length_sq,          length,             length_sq_per_time>,
    std::tuple<length_sq_per_time, length_sq_per_time, length_sq_per_time_sq>
>;

int test_main (int, char *[])
{
    static_assert(bub::detail::has_identity<E_matrix_2x2_double_type>::value, "");
    static_assert(bub::detail::has_identity<derived_from_E_matrix_2x2_double_type>::value, "");
    static_assert(bub::detail::has_identity<E_matrix_2x2_time_type>::value, "");
    static_assert(bub::detail::has_identity<derived_from_E_matrix_2x2_time_type>::value, "");

    static_assert(!bub::detail::has_identity<int>::value, "");
    static_assert(!bub::detail::has_identity<uninvertible_matrix_2x2>::value, "");
    static_assert(!bub::detail::has_identity<uninvertible_matrix_3x3>::value, "");

    return 0;
}
