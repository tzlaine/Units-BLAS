// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/result_of/sum.hpp>

#include <boost/test/minimal.hpp>


int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::sum<A_matrix_2x1_float_int_type>::type,
            float
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::sum<C_matrix_1x2_int_type>::type,
            int
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::sum<D_matrix_1x2_int_double_type>::type,
            double
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::sum<B_matrix_2x1_length_type>::type,
            length
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::sum<C_matrix_1x2_length_type>::type,
            length
        >::type
    ));

    return 0;
}
