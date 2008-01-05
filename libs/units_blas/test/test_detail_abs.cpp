// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/detail/abs.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/length.hpp>

#include <boost/test/minimal.hpp>


typedef boost::units::quantity<boost::units::SI::length, float> length;

namespace bub = boost::units_blas;

int test_main (int, char *[])
{
    int neg_i = -1;
    int pos_i = 1;

    float neg_f = -1.0f;
    float pos_f = 1.0f;

    double neg_d = -1.0;
    double pos_d = 1.0;

    length neg_lqt = length::from_value(-1.0f);
    length pos_lqt = length::from_value(1.0f);

    BOOST_CHECK((bub::detail::abs_(neg_i) == pos_i));
    BOOST_CHECK((bub::detail::abs_(pos_i) == pos_i));
    BOOST_CHECK((bub::detail::abs_(neg_f) == pos_f));
    BOOST_CHECK((bub::detail::abs_(pos_f) == pos_f));
    BOOST_CHECK((bub::detail::abs_(neg_d) == pos_d));
    BOOST_CHECK((bub::detail::abs_(pos_d) == pos_d));
    BOOST_CHECK((bub::detail::abs_(neg_lqt) == pos_lqt));
    BOOST_CHECK((bub::detail::abs_(pos_lqt) == pos_lqt));

    return 0;
}
