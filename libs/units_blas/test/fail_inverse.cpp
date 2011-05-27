// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/operations.hpp>


typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<length_sq,          length_sq,          length_sq_per_time>,
        boost::fusion::vector<length_sq,          length_sq,          length_sq_per_time>,
        boost::fusion::vector<length_sq_per_time, length_sq_per_time, length_sq_per_time_sq>
    >
>::type matrix_3x3_units_type;

int test_main (int, char *[])
{
    inverse(matrix_3x3_units_type());

    return 0;
}
