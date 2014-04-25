// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/make_matrix.hpp>
#include <boost/units_blas/io.hpp>

#include <boost/units/io.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/length.hpp>

#include <boost/test/minimal.hpp>

#include <iostream>


typedef boost::units::quantity<boost::units::si::length, float> length;

namespace bub = boost::units_blas;

typedef bub::uniform_matrix<
    float, 4, 1
> matrix_4x1_float_type;

typedef bub::uniform_matrix<
    float, 1, 4
> matrix_1x4_float_type;

typedef bub::uniform_matrix<
    float, 4, 4
> matrix_4x4_float_type;

typedef bub::uniform_matrix<
    length, 4, 1
> matrix_4x1_length_type;

typedef bub::uniform_matrix<
    length, 1, 4
> matrix_1x4_length_type;

typedef bub::uniform_matrix<
    length, 4, 4
> matrix_4x4_length_type;

int test_main (int, char *[])
{
    std::cout
        << matrix_4x1_float_type{}
        << matrix_1x4_float_type{}
        << matrix_4x4_float_type{}
        << matrix_4x1_length_type{}
        << matrix_1x4_length_type{}
        << matrix_4x4_length_type{}
        ;

    return 0;
}
