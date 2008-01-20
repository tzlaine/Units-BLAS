// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "failure_tests.hpp"

#include <boost/units_blas/operations.hpp>

#include <boost/test/minimal.hpp>


int test_main (int, char *[])
{
    matrix_2x2 m1;
    matrix_2x1 m2;
    m2 = m1;

    return 0;
}
