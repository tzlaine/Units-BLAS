// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/result_of/cross_product.hpp>

#include <boost/test/minimal.hpp>


typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int>,
        boost::fusion::vector<float>,
        boost::fusion::vector<double>
    >
> F_matrix_3x1_mixed_fundamentals_type;

typedef bub::result_of::cross_product<
    F_matrix_3x1_mixed_fundamentals_type,
    F_matrix_3x1_mixed_fundamentals_type
>::type F_cross_F;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<double>,
        boost::fusion::vector<double>,
        boost::fusion::vector<float>
    >
> F_cross_F_manual;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, int, long>
    >
> G_matrix_1x3_mixed_fundamentals_type_1;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, int, int>
    >
> H_matrix_1x3_mixed_fundamentals_type_2;

typedef bub::result_of::cross_product<
    G_matrix_1x3_mixed_fundamentals_type_1,
    H_matrix_1x3_mixed_fundamentals_type_2
>::type G_cross_H;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<long, long, int>
    >
> G_cross_H_manual;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length>,
        boost::fusion::vector<length>,
        boost::fusion::vector<length>
    >
> I_matrix_3x1_mixed_units_type;

typedef bub::result_of::cross_product<
    I_matrix_3x1_mixed_units_type,
    I_matrix_3x1_mixed_units_type
>::type I_cross_I;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length_sq>,
        boost::fusion::vector<length_sq>,
        boost::fusion::vector<length_sq>
    >
> I_cross_I_manual;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_, length, frequency>
    >
> J_matrix_1x3_mixed_units_type;

typedef bub::result_of::cross_product<
    J_matrix_1x3_mixed_units_type,
    J_matrix_1x3_mixed_units_type
>::type J_cross_J;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<velocity, dimensionless, time_length>
    >
> J_cross_J_manual;

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((boost::is_same<F_cross_F, F_cross_F_manual>));
    BOOST_MPL_ASSERT((boost::is_same<G_cross_H, G_cross_H_manual>));
    BOOST_MPL_ASSERT((boost::is_same<I_cross_I, I_cross_I_manual>));
    BOOST_MPL_ASSERT((boost::is_same<J_cross_J, J_cross_J_manual>));

    return 0;
}
