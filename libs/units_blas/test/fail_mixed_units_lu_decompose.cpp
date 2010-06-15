// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/operations.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/cgs/length.hpp>

#include <boost/test/minimal.hpp>


namespace bub = boost::units_blas;

typedef boost::units::quantity<boost::units::si::length> length;
typedef boost::units::quantity<boost::units::cgs::length> length2;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, length2>,
        boost::fusion::vector<length, length>
    >
> A;

int test_main (int, char *[])
{
    bub::inverse(A());

    return 0;
}
