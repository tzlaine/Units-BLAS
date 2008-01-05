// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/result_of/transpose.hpp>

#include <boost/test/minimal.hpp>


typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<float, int>
    >
> A_matrix_2x1_float_int_type_transpose;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int>,
        boost::fusion::vector<double>
    >
> D_matrix_1x2_int_double_type_transpose;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_, length>
    >
> A_matrix_2x1_time_length_type_transpose;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length>,
        boost::fusion::vector<time_>
    >
> D_matrix_1x2_length_time_type_transpose;

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<A_matrix_2x1_float_int_type>::type,
            A_matrix_2x1_float_int_type_transpose
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<B_matrix_2x1_int_type>::type,
            C_matrix_1x2_int_type
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<C_matrix_1x2_int_type>::type,
            B_matrix_2x1_int_type
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<D_matrix_1x2_int_double_type>::type,
            D_matrix_1x2_int_double_type_transpose
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<E_matrix_2x2_double_type>::type,
            E_matrix_2x2_double_type
        >
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<A_matrix_2x1_time_length_type>::type,
            A_matrix_2x1_time_length_type_transpose
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<B_matrix_2x1_length_type>::type,
            C_matrix_1x2_length_type
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<C_matrix_1x2_length_type>::type,
            B_matrix_2x1_length_type
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<D_matrix_1x2_length_time_type>::type,
            D_matrix_1x2_length_time_type_transpose
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<E_matrix_2x2_time_type>::type,
            E_matrix_2x2_time_type
        >
    ));


    // round-trip test

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::transpose<bub::result_of::transpose<A_matrix_2x1_float_int_type>::type>::type,
            A_matrix_2x1_float_int_type
        >
    ));

    return 0;
}
