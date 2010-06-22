// Copyright (C) 2010 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units_blas/dont_care.hpp>

#include <boost/test/minimal.hpp>


namespace bub = boost::units_blas;

typedef boost::units::quantity<boost::units::si::length> length;

int test_main (int, char *[])
{
    int i;
    double d;
    length l;

    bub::_ _;

    bub::_ _2(i);
    bub::_ _3(d);
    bub::_ _4(l);

    static_cast<bub::_>(i);
    static_cast<bub::_>(d);
    static_cast<bub::_>(l);

    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ * _)), bub::_>::type));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((i * _)), bub::_>::type));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ * i)), bub::_>::type));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((d * _)), bub::_>::type));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ * d)), bub::_>::type));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((l * _)), bub::_>::type));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ * l)), bub::_>::type));

    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ + _)), bub::_>::type));

    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((i + _)), int>::type));
    BOOST_CHECK((i + _ == i));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ + i)), int>::type));
    BOOST_CHECK((_ + i == i));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((d + _)), double>::type));
    BOOST_CHECK((d + _ == d));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ + d)), double>::type));
    BOOST_CHECK((_ + d == d));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((l + _)), length>::type));
    BOOST_CHECK((l + _ == l));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ + l)), length>::type));
    BOOST_CHECK((_ + l == l));

    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ - _)), bub::_>::type));

    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((i - _)), int>::type));
    BOOST_CHECK((i - _ == i));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ - i)), int>::type));
    BOOST_CHECK((_ - i == i));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((d - _)), double>::type));
    BOOST_CHECK((d - _ == d));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ - d)), double>::type));
    BOOST_CHECK((_ - d == d));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((l - _)), length>::type));
    BOOST_CHECK((l - _ == l));
    BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((_ - l)), length>::type));
    BOOST_CHECK((_ - l == l));

    return 0;
}
