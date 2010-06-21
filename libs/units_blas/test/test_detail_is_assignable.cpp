// Copyright (C) 2010 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/detail/is_assignable.hpp>

#include <boost/test/minimal.hpp>


typedef boost::units::quantity<boost::units::si::length, float> float_length;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<float_length>,
        boost::fusion::vector<float_length>
    >
> matrix_2x1_float_length_type;

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((bub::detail::is_assignable<A_matrix_2x1_float_int_type, A_matrix_2x1_float_int_type>));

    BOOST_MPL_ASSERT((bub::detail::is_assignable<A_matrix_2x1_float_int_type, B_matrix_2x1_int_type>));
    BOOST_MPL_ASSERT((bub::detail::is_assignable<B_matrix_2x1_int_type, A_matrix_2x1_float_int_type>));

    BOOST_MPL_ASSERT_NOT((bub::detail::is_assignable<B_matrix_2x1_length_type, B_matrix_2x1_int_type>));
    BOOST_MPL_ASSERT_NOT((bub::detail::is_assignable<B_matrix_2x1_int_type, B_matrix_2x1_length_type>));

    BOOST_MPL_ASSERT((bub::detail::is_assignable<matrix_2x1_float_length_type, B_matrix_2x1_length_type>));

    // TODO: This should work, but currently does not.
    //BOOST_MPL_ASSERT((bub::detail::is_assignable<B_matrix_2x1_length_type, matrix_2x1_float_length_type>));

    BOOST_MPL_ASSERT_NOT((bub::detail::is_assignable<B_matrix_2x1_int_type, C_matrix_1x2_int_type>));
    BOOST_MPL_ASSERT_NOT((bub::detail::is_assignable<C_matrix_1x2_int_type, B_matrix_2x1_int_type>));

    return 0;
}
