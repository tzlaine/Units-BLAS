// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/result_of/dot_product.hpp>

#include <boost/test/minimal.hpp>


int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::dot_product<
                A_matrix_2x1_float_int_type,
                B_matrix_2x1_int_type
            >::type,
            float
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::dot_product<
                B_matrix_2x1_int_type,
                A_matrix_2x1_float_int_type
            >::type,
            float
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::dot_product<
                bub::result_of::transpose<C_matrix_1x2_int_type>::type,
                bub::result_of::transpose<D_matrix_1x2_int_double_type>::type
            >::type,
            double
        >
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::dot_product<
                bub::result_of::transpose<D_matrix_1x2_int_double_type>::type,
                bub::result_of::transpose<C_matrix_1x2_int_type>::type
            >::type,
            double
        >
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::dot_product<
                B_matrix_2x1_length_type,
                bub::result_of::transpose<C_matrix_1x2_length_type>::type
            >::type,
            length_sq
        >
    ));

    return 0;
}
