// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/detail/zero_value.hpp>
#include <boost/units_blas/detail/one_value.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/length.hpp>

#include <boost/test/minimal.hpp>


typedef boost::units::quantity<boost::units::si::length, float> length;

namespace bub = boost::units_blas;

int test_main (int, char *[])
{
    BOOST_CHECK((bub::detail::zero_value<int>::value() == 0));
    BOOST_CHECK((bub::detail::zero_value<float>::value() == 0.0f));
    BOOST_CHECK((bub::detail::zero_value<double>::value() == 0.0));
    BOOST_CHECK((bub::detail::zero_value<length>::value() == length::from_value(0.0f)));

    BOOST_CHECK((bub::detail::one_value<int>::value() == 1));
    BOOST_CHECK((bub::detail::one_value<float>::value() == 1.0f));
    BOOST_CHECK((bub::detail::one_value<double>::value() == 1.0));
    BOOST_CHECK((bub::detail::one_value<length>::value() == length::from_value(1.0f)));

    return 0;
}
