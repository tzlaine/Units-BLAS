// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/make_matrix.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/length.hpp>

#include <boost/test/minimal.hpp>


namespace bub = boost::units_blas;
namespace bh = boost::hana;

typedef bub::matrix<
    bh::_tuple<float>,
    bh::_tuple<float>
> manual_matrix_2x1_type;

typedef bub::matrix<
    bh::_tuple<float, float>
> manual_matrix_1x2_type;

typedef bub::matrix<
    bh::_tuple<float, float>,
    bh::_tuple<float, float>
> manual_matrix_2x2_type;

typedef boost::units::quantity<boost::units::si::length> length;

typedef bub::matrix<
    bh::_tuple<length>,
    bh::_tuple<length>
> manual_matrix_2x1_units_type;

typedef bub::matrix<
    bh::_tuple<length, length>
> manual_matrix_1x2_units_type;

typedef bub::matrix<
    bh::_tuple<length, length>,
    bh::_tuple<length, length>
> manual_matrix_2x2_units_type;

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<float, 1, 2>, manual_matrix_1x2_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<float, 2, 1>, manual_matrix_2x1_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<float, 2, 2>, manual_matrix_2x2_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<length, 1, 2>, manual_matrix_1x2_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<length, 2, 1>, manual_matrix_2x1_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<length, 2, 2>, manual_matrix_2x2_units_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::vector<float, float>, manual_matrix_2x1_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_vector<float, 2>, manual_matrix_2x1_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::vector<length, length>, manual_matrix_2x1_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_vector<length, 2>, manual_matrix_2x1_units_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<float, float>, manual_matrix_1x2_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_transpose_vector<float, 2>, manual_matrix_1x2_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<length, length>, manual_matrix_1x2_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_transpose_vector<length, 2>, manual_matrix_1x2_units_type>));

    return 0;
}
