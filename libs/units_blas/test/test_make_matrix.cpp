// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/make_matrix.hpp>

#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/list.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/list.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/length.hpp>

#include <boost/test/minimal.hpp>


namespace bub = boost::units_blas;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<float>,
        boost::fusion::vector<float>
    >
> manual_matrix_2x1_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<float, float>
    >
> manual_matrix_1x2_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<float, float>,
        boost::fusion::vector<float, float>
    >
> manual_matrix_2x2_type;

typedef bub::make_matrix<
    boost::mpl::vector<
        boost::mpl::vector<float>,
        boost::mpl::vector<float>
    >
>::type made_matrix_2x1_type;

typedef bub::make_matrix<
    boost::mpl::vector<
        boost::mpl::vector<float, float>
    >
>::type made_matrix_1x2_type;

typedef bub::make_matrix<
    boost::mpl::vector<
        boost::mpl::vector<float, float>,
        boost::mpl::vector<float, float>
    >
>::type made_matrix_2x2_type;

typedef boost::units::quantity<boost::units::SI::length> length;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length>,
        boost::fusion::vector<length>
    >
> manual_matrix_2x1_units_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, length>
    >
> manual_matrix_1x2_units_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, length>,
        boost::fusion::vector<length, length>
    >
> manual_matrix_2x2_units_type;

typedef bub::make_matrix<
    boost::mpl::vector<
        boost::mpl::vector<length>,
        boost::mpl::vector<length>
    >
>::type made_matrix_2x1_units_type;

typedef bub::make_matrix<
    boost::mpl::vector<
        boost::mpl::vector<length, length>
    >
>::type made_matrix_1x2_units_type;

typedef bub::make_matrix<
    boost::mpl::vector<
        boost::mpl::vector<length, length>,
        boost::mpl::vector<length, length>
    >
>::type made_matrix_2x2_units_type;

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((boost::is_same<manual_matrix_1x2_type, made_matrix_1x2_type>));
    BOOST_MPL_ASSERT((boost::is_same<manual_matrix_2x1_type, made_matrix_2x1_type>));
    BOOST_MPL_ASSERT((boost::is_same<manual_matrix_2x2_type, made_matrix_2x2_type>));

    BOOST_MPL_ASSERT((boost::is_same<manual_matrix_1x2_units_type, made_matrix_1x2_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<manual_matrix_2x1_units_type, made_matrix_2x1_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<manual_matrix_2x2_units_type, made_matrix_2x2_units_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<float, 1, 2>::type, manual_matrix_1x2_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<float, 2, 1>::type, manual_matrix_2x1_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<float, 2, 2>::type, manual_matrix_2x2_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<length, 1, 2>::type, manual_matrix_1x2_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<length, 2, 1>::type, manual_matrix_2x1_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_matrix<length, 2, 2>::type, manual_matrix_2x2_units_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::vector<boost::fusion::vector<float, float> >::type, manual_matrix_2x1_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::vector<boost::fusion::list<float, float> >::type, manual_matrix_2x1_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::vector<boost::mpl::vector<float, float> >::type, manual_matrix_2x1_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::vector<boost::mpl::list<float, float> >::type, manual_matrix_2x1_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_vector<float, 2>::type, manual_matrix_2x1_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::vector<boost::fusion::vector<length, length> >::type, manual_matrix_2x1_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::vector<boost::fusion::list<length, length> >::type, manual_matrix_2x1_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::vector<boost::mpl::vector<length, length> >::type, manual_matrix_2x1_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::vector<boost::mpl::list<length, length> >::type, manual_matrix_2x1_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_vector<length, 2>::type, manual_matrix_2x1_units_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<boost::fusion::vector<float, float> >::type, manual_matrix_1x2_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<boost::fusion::list<float, float> >::type, manual_matrix_1x2_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<boost::mpl::vector<float, float> >::type, manual_matrix_1x2_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<boost::mpl::list<float, float> >::type, manual_matrix_1x2_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_transpose_vector<float, 2>::type, manual_matrix_1x2_type>));

    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<boost::fusion::vector<length, length> >::type, manual_matrix_1x2_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<boost::fusion::list<length, length> >::type, manual_matrix_1x2_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<boost::mpl::vector<length, length> >::type, manual_matrix_1x2_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::transpose_vector<boost::mpl::list<length, length> >::type, manual_matrix_1x2_units_type>));
    BOOST_MPL_ASSERT((boost::is_same<bub::uniform_transpose_vector<length, 2>::type, manual_matrix_1x2_units_type>));

    return 0;
}
