// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/result_of/inverse.hpp>

#include <boost/test/minimal.hpp>


typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, long>,
        boost::fusion::vector<float, double>
    >
> matrix_2x2_mixed_fundamentals_type;

typedef bub::result_of::inverse<
    matrix_2x2_mixed_fundamentals_type
>::type matrix_2x2_mixed_fundamentals_type_inverse;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_, frequency>,
        boost::fusion::vector<time_, frequency>
    >
> matrix_2x2_mixed_units_type;

typedef bub::result_of::inverse<
    matrix_2x2_mixed_units_type
>::type matrix_2x2_mixed_units_type_inverse;

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_2x2_mixed_fundamentals_type_inverse, 0, 0>::type,
            int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_2x2_mixed_fundamentals_type_inverse, 0, 1>::type,
            float
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_2x2_mixed_fundamentals_type_inverse, 1, 0>::type,
            long
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_2x2_mixed_fundamentals_type_inverse, 1, 1>::type,
            double
        >::type
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_2x2_mixed_units_type_inverse, 0, 0>::type,
            frequency
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_2x2_mixed_units_type_inverse, 0, 1>::type,
            frequency
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_2x2_mixed_units_type_inverse, 1, 0>::type,
            time_
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_2x2_mixed_units_type_inverse, 1, 1>::type,
            time_
        >::type
    ));

    return 0;
}
