// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/result_of/determinant.hpp>

#include <boost/test/minimal.hpp>


typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int>
    >
> matrix_1x1_fundamentals_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, int>,
        boost::fusion::vector<int, float>
    >
> matrix_2x2_fundamentals_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, int, int>,
        boost::fusion::vector<int, double, float>,
        boost::fusion::vector<int, int, int>
    >
> matrix_3x3_fundamentals_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, int, int, int>,
        boost::fusion::vector<int, int, float, int>,
        boost::fusion::vector<int, int, int, int>,
        boost::fusion::vector<int, int, int, int>
    >
> matrix_4x4_fundamentals_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length>
    >
> matrix_1x1_units_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, length>,
        boost::fusion::vector<time_, time_>
    >
> matrix_2x2_units_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, length, length>,
        boost::fusion::vector<frequency, frequency, frequency>,
        boost::fusion::vector<time_, time_, time_>
    >
> matrix_3x3_units_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<frequency, frequency, frequency, frequency>,
        boost::fusion::vector<length, length, length, length>,
        boost::fusion::vector<frequency, frequency, frequency, frequency>,
        boost::fusion::vector<time_, time_, time_, time_>
    >
> matrix_4x4_units_type;

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::determinant<matrix_1x1_fundamentals_type>::type,
            int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::determinant<matrix_2x2_fundamentals_type>::type,
            float
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::determinant<matrix_3x3_fundamentals_type>::type,
            double
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::determinant<matrix_4x4_fundamentals_type>::type,
            float
        >::type
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::determinant<matrix_1x1_units_type>::type,
            length
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::determinant<matrix_2x2_units_type>::type,
            time_length
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::determinant<matrix_3x3_units_type>::type,
            length
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::determinant<matrix_4x4_units_type>::type,
            velocity
        >::type
    ));

    return 0;
}
