// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/operations.hpp>


using uninvertible_matrix_3x3 = bub::matrix<
    std::tuple<length_sq,          length_sq,          length_sq_per_time>,
    std::tuple<length_sq,          length,             length_sq_per_time>,
    std::tuple<length_sq_per_time, length_sq_per_time, length_sq_per_time_sq>
>;

int test_main (int, char *[])
{
    inverse(uninvertible_matrix_3x3{});

    return 0;
}
