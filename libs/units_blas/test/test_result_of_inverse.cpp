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
        boost::fusion::vector<float>
    >
> matrix_1x1_fundamentals_type;

typedef bub::result_of::inverse<
    matrix_1x1_fundamentals_type
>::type matrix_1x1_fundamentals_type_inverse;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_>
    >
> matrix_1x1_units_type;

typedef bub::result_of::inverse<
    matrix_1x1_units_type
>::type matrix_1x1_units_type_inverse;

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

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, long, char>,
        boost::fusion::vector<short, unsigned int, unsigned long>,
        boost::fusion::vector<float, double, unsigned char>
    >
> matrix_3x3_mixed_fundamentals_type;

typedef bub::result_of::inverse<
    matrix_3x3_mixed_fundamentals_type
>::type matrix_3x3_mixed_fundamentals_type_inverse;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length_sq, length_sq, length_sq_per_time>,
        boost::fusion::vector<length_sq, length_sq, length_sq_per_time>,
        boost::fusion::vector<length_sq_per_time, length_sq_per_time, length_sq_per_time_sq>
    >
> matrix_3x3_mixed_units_type;

typedef bub::result_of::inverse<
    matrix_3x3_mixed_units_type
>::type matrix_3x3_mixed_units_type_inverse;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, long, char, int>,
        boost::fusion::vector<short, unsigned int, unsigned long, int>,
        boost::fusion::vector<float, double, unsigned char, int>,
        boost::fusion::vector<float, double, int, int>
    >
> matrix_4x4_mixed_fundamentals_type;

typedef bub::result_of::inverse<
    matrix_4x4_mixed_fundamentals_type
>::type matrix_4x4_mixed_fundamentals_type_inverse;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length_sq, length_sq, length_sq, length_sq_per_time>,
        boost::fusion::vector<length_sq, length_sq, length_sq, length_sq_per_time>,
        boost::fusion::vector<length_sq, length_sq, length_sq, length_sq_per_time>,
        boost::fusion::vector<length_sq_per_time, length_sq_per_time, length_sq_per_time, length_sq_per_time_sq>
    >
> matrix_4x4_mixed_units_type;

typedef bub::result_of::inverse<
    matrix_4x4_mixed_units_type
>::type matrix_4x4_mixed_units_type_inverse;

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_1x1_fundamentals_type_inverse, 0, 0>::type,
            float
        >::type
    ));


    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_1x1_units_type_inverse, 0, 0>::type,
            frequency
        >::type
    ));


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


    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_fundamentals_type_inverse, 0, 0>::type,
            int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_fundamentals_type_inverse, 0, 1>::type,
            short
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_fundamentals_type_inverse, 0, 2>::type,
            float
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_fundamentals_type_inverse, 1, 0>::type,
            long
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_fundamentals_type_inverse, 1, 1>::type,
            unsigned int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_fundamentals_type_inverse, 1, 2>::type,
            double
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_fundamentals_type_inverse, 2, 0>::type,
            char
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_fundamentals_type_inverse, 2, 1>::type,
            unsigned long
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_fundamentals_type_inverse, 2, 2>::type,
            unsigned char
        >::type
    ));


    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_units_type_inverse, 0, 0>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_units_type_inverse, 0, 1>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_units_type_inverse, 0, 2>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_units_type_inverse, 1, 0>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_units_type_inverse, 1, 1>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_units_type_inverse, 1, 2>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_units_type_inverse, 2, 0>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_units_type_inverse, 2, 1>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_3x3_mixed_units_type_inverse, 2, 2>::type,
            time_sq_per_length_sq
        >::type
    ));


    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 0, 0>::type,
            int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 0, 1>::type,
            short
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 0, 2>::type,
            float
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 0, 3>::type,
            float
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 1, 0>::type,
            long
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 1, 1>::type,
            unsigned int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 1, 2>::type,
            double
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 1, 3>::type,
            double
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 2, 0>::type,
            char
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 2, 1>::type,
            unsigned long
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 2, 2>::type,
            unsigned char
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 2, 3>::type,
            int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 3, 0>::type,
            int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 3, 1>::type,
            int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 3, 2>::type,
            int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_fundamentals_type_inverse, 3, 3>::type,
            int
        >::type
    ));


    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 0, 0>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 0, 1>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 0, 2>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 0, 3>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 1, 0>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 1, 1>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 1, 2>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 1, 3>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 2, 0>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 2, 1>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 2, 2>::type,
            inv_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 2, 3>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 3, 0>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 3, 1>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 3, 2>::type,
            time_per_length_sq
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::value_at_c<matrix_4x4_mixed_units_type_inverse, 3, 3>::type,
            time_sq_per_length_sq
        >::type
    ));

    return 0;
}
