// Copyright (C) 2010 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/dont_care.hpp>
#include <boost/units_blas/matrix.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/length.hpp>


namespace bub = boost::units_blas;

int test_main (int, char *[])
{
    typedef boost::units::quantity<boost::units::si::time> time_;
    typedef boost::units::quantity<boost::units::si::length> length;

    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_>,
            boost::fusion::vector<bub::_>
        >
    > matrix_2x1_time_dont_care;

    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_>,
            boost::fusion::vector<length>
        >
    > matrix_2x1_time_length;

    matrix_2x1_time_length time_length = matrix_2x1_time_dont_care();

    return 0;
}
