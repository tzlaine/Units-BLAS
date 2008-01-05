// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/detail/get_value_type.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/length.hpp>

#include <boost/test/minimal.hpp>


struct value_type_haver
{
    typedef int value_type;
};

typedef boost::units::quantity<boost::units::SI::length, float> length;

namespace bub = boost::units_blas;

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((bub::detail::has_value_type<value_type_haver>));
    BOOST_MPL_ASSERT((bub::detail::has_value_type<length>));
    BOOST_MPL_ASSERT_NOT((bub::detail::has_value_type<int>));

    BOOST_MPL_ASSERT((boost::mpl::is_same<bub::detail::get_value_type<int>::type, int>));
    BOOST_MPL_ASSERT((boost::mpl::is_same<bub::detail::get_value_type<value_type_haver>::type, int>));
    BOOST_MPL_ASSERT((boost::mpl::is_same<bub::detail::get_value_type<length>::type, float>));

    return 0;
}
