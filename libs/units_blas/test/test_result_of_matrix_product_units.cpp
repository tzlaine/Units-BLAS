// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/traits.hpp>
#include <boost/units_blas/result_of/matrix_product.hpp>
#include <boost/units_blas/result_of/value_at.hpp>

#include <boost/test/minimal.hpp>


int test_main (int, char *[])
{
    // matrix-matrix products

    typedef bub::result_of::matrix_product<
        A_matrix_2x1_time_length_type,
        C_matrix_1x2_length_type
    >::type AxC;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_length, time_length>,
            boost::fusion::vector<length_sq, length_sq>
        >
    >::type AxC_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<AxC, AxC_manual_product>));

    typedef bub::result_of::matrix_product<
        A_matrix_2x1_time_length_type,
        D_matrix_1x2_length_time_type
    >::type AxD;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_length, time_sq>,
            boost::fusion::vector<length_sq, time_length>
        >
    >::type AxD_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<AxD, AxD_manual_product>));

    typedef bub::result_of::matrix_product<
        D_matrix_1x2_length_time_type,
        A_matrix_2x1_time_length_type
    >::type DxA;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_length>
        >
    >::type DxA_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<DxA, DxA_manual_product>));

    typedef bub::result_of::matrix_product<
        B_matrix_2x1_length_type,
        C_matrix_1x2_length_type
    >::type BxC;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<length_sq, length_sq>,
            boost::fusion::vector<length_sq, length_sq>
        >
    >::type BxC_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<BxC, BxC_manual_product>));

    typedef bub::result_of::matrix_product<
        C_matrix_1x2_length_type,
        B_matrix_2x1_length_type
    >::type CxB;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<length_sq>
        >
    >::type CxB_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<CxB, CxB_manual_product>));

    typedef bub::result_of::matrix_product<
        B_matrix_2x1_length_type,
        D_matrix_1x2_length_time_type
    >::type BxD;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<length_sq, time_length>,
            boost::fusion::vector<length_sq, time_length>
        >
    >::type BxD_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<BxD, BxD_manual_product>));

    typedef bub::result_of::matrix_product<
        C_matrix_1x2_length_type,
        E_matrix_2x2_time_type
    >::type CxE;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_length, time_length>
        >
    >::type CxE_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<CxE, CxE_manual_product>));

    typedef bub::result_of::matrix_product<
        E_matrix_2x2_time_type,
        B_matrix_2x1_length_type
    >::type ExB;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_length>,
            boost::fusion::vector<time_length>
        >
    >::type ExB_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<ExB, ExB_manual_product>));

    typedef bub::result_of::matrix_product<
        E_matrix_2x2_time_type,
        E_matrix_2x2_time_type
    >::type ExE;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_sq, time_sq>,
            boost::fusion::vector<time_sq, time_sq>
        >
    >::type ExE_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<ExE, ExE_manual_product>));


    // matrix-scalar products

    typedef bub::result_of::scalar_product<
        A_matrix_2x1_time_length_type,
        double
    >::type Axdouble;
    BOOST_MPL_ASSERT((boost::is_same<Axdouble, A_matrix_2x1_time_length_type>));

    typedef bub::result_of::scalar_product<
        A_matrix_2x1_time_length_type,
        time_
    >::type Axtime;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_sq>,
            boost::fusion::vector<time_length>
        >
    >::type Axtime_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<Axtime, Axtime_manual_product>));

    typedef bub::result_of::scalar_product<
        D_matrix_1x2_length_time_type,
        double
    >::type Dxdouble;
    BOOST_MPL_ASSERT((boost::is_same<Dxdouble, D_matrix_1x2_length_time_type>));

    typedef bub::result_of::scalar_product<
        D_matrix_1x2_length_time_type,
        time_
    >::type Dxtime;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_length, time_sq>
        >
    >::type Dxtime_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<Dxtime, Dxtime_manual_product>));


    // matrix-scalar quotients

    typedef bub::result_of::scalar_quotient<
        A_matrix_2x1_time_length_type,
        double
    >::type Addouble;
    BOOST_MPL_ASSERT((boost::is_same<Addouble, A_matrix_2x1_time_length_type>));

    typedef bub::result_of::scalar_quotient<
        A_matrix_2x1_time_length_type,
        time_
    >::type Adtime;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<dimensionless>,
            boost::fusion::vector<velocity>
        >
    >::type Adtime_manual_quotient;
    BOOST_MPL_ASSERT((boost::is_same<Adtime, Adtime_manual_quotient>));

    typedef bub::result_of::scalar_quotient<
        D_matrix_1x2_length_time_type,
        double
    >::type Dddouble;
    BOOST_MPL_ASSERT((boost::is_same<Dddouble, D_matrix_1x2_length_time_type>));

    typedef bub::result_of::scalar_quotient<
        D_matrix_1x2_length_time_type,
        time_
    >::type Ddtime;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<velocity, dimensionless>
        >
    >::type Ddtime_manual_quotient;
    BOOST_MPL_ASSERT((boost::is_same<Ddtime, Ddtime_manual_quotient>));


    // element-wise matrix sums

    typedef boost::units::quantity<boost::units::si::time, float> time_float;
    typedef boost::units::quantity<boost::units::si::length, float> length_float;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_float>,
            boost::fusion::vector<length_float>
       >
    >::type A_matrix_2x1_time_length_float_type;

    typedef bub::result_of::matrix_element_sum<
        A_matrix_2x1_time_length_type,
        A_matrix_2x1_time_length_float_type
    >::type ApB_elementwise_sum;
    BOOST_MPL_ASSERT((boost::is_same<ApB_elementwise_sum, A_matrix_2x1_time_length_type>));


    // element-wise matrix differences

    typedef bub::result_of::matrix_element_difference<
        A_matrix_2x1_time_length_type,
        A_matrix_2x1_time_length_float_type
    >::type AmB_elementwise_difference;
    BOOST_MPL_ASSERT((boost::is_same<AmB_elementwise_difference, A_matrix_2x1_time_length_type>));


    // element-wise matrix products

    typedef bub::result_of::matrix_element_product<
        A_matrix_2x1_time_length_type,
        B_matrix_2x1_length_type
    >::type AxB_elementwise_product;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_length>,
            boost::fusion::vector<length_sq>
        >
    >::type AxB_elementwise_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<AxB_elementwise_product, AxB_elementwise_manual_product>));

    typedef bub::result_of::matrix_element_product<
        C_matrix_1x2_length_type,
        D_matrix_1x2_length_time_type
    >::type CxD_elementwise_product;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<length_sq, time_length>
        >
    >::type CxD_elementwise_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<CxD_elementwise_product, CxD_elementwise_manual_product>));


    // element-wise matrix quotients

    typedef bub::result_of::matrix_element_quotient<
        B_matrix_2x1_length_type,
        A_matrix_2x1_time_length_type
    >::type BxA_elementwise_quotient;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<velocity>,
            boost::fusion::vector<dimensionless>
        >
    >::type BxA_elementwise_manual_quotient;
    BOOST_MPL_ASSERT((boost::is_same<BxA_elementwise_quotient, BxA_elementwise_manual_quotient>));

    typedef bub::result_of::matrix_element_quotient<
        C_matrix_1x2_length_type,
        D_matrix_1x2_length_time_type
    >::type CxD_elementwise_quotient;
    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<dimensionless, velocity>
        >
    >::type CxD_elementwise_manual_quotient;
    BOOST_MPL_ASSERT((boost::is_same<CxD_elementwise_quotient, CxD_elementwise_manual_quotient>));

    return 0;
}
