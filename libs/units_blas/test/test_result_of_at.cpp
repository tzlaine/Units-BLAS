// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/result_of/at.hpp>

#include <boost/test/minimal.hpp>


int test_main (int, char *[])
{
    typedef bub::size_t_<0> zero;
    typedef bub::size_t_<1> one;

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<A_matrix_2x1_float_int_type, zero, zero>::type,
            float &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<A_matrix_2x1_float_int_type, one, zero>::type,
            int &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<A_matrix_2x1_float_int_type, 0, 0>::type,
            float &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<A_matrix_2x1_float_int_type, 1, 0>::type,
            int &
        >
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<A_matrix_2x1_float_int_type const, one, zero>::type,
            int const &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<A_matrix_2x1_float_int_type const, 1, 0>::type,
            int const &
        >
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<A_matrix_2x1_time_length_type, zero, zero>::type,
            time_ &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<A_matrix_2x1_time_length_type, one, zero>::type,
            length &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<A_matrix_2x1_time_length_type, 0, 0>::type,
            time_ &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<A_matrix_2x1_time_length_type, 1, 0>::type,
            length &
        >
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<A_matrix_2x1_time_length_type const, one, zero>::type,
            length const &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<A_matrix_2x1_time_length_type const, 1, 0>::type,
            length const &
        >
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<D_matrix_1x2_int_double_type, zero, zero>::type,
            int &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<D_matrix_1x2_int_double_type, zero, one>::type,
            double &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<D_matrix_1x2_int_double_type, 0, 0>::type,
            int &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<D_matrix_1x2_int_double_type, 0, 1>::type,
            double &
        >
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<D_matrix_1x2_length_time_type, zero, zero>::type,
            length &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at<D_matrix_1x2_length_time_type, zero, one>::type,
            time_ &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<D_matrix_1x2_length_time_type, 0, 0>::type,
            length &
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::at_c<D_matrix_1x2_length_time_type, 0, 1>::type,
            time_ &
        >
    ));

    return 0;
}
