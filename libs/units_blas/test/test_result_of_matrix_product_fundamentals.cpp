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
        A_matrix_2x1_float_int_type,
        C_matrix_1x2_int_type
    >::type AxC;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<float, float>,
            boost::fusion::vector<int, int>
        >
    > AxC_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<AxC, AxC_manual_product>));

    typedef bub::result_of::matrix_product<
        C_matrix_1x2_int_type,
        A_matrix_2x1_float_int_type
    >::type CxA;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<float>
        >
    > CxA_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<CxA, CxA_manual_product>));

    typedef bub::result_of::matrix_product<
        A_matrix_2x1_float_int_type,
        D_matrix_1x2_int_double_type
    >::type AxD;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<float, double>,
            boost::fusion::vector<int, double>
        >
    > AxD_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<AxD, AxD_manual_product>));

    typedef bub::result_of::matrix_product<
        D_matrix_1x2_int_double_type,
        A_matrix_2x1_float_int_type
    >::type DxA;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<double>
        >
    > DxA_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<DxA, DxA_manual_product>));

    typedef bub::result_of::matrix_product<
        B_matrix_2x1_int_type,
        C_matrix_1x2_int_type
    >::type BxC;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<int, int>,
            boost::fusion::vector<int, int>
        >
    > BxC_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<BxC, BxC_manual_product>));

    typedef bub::result_of::matrix_product<
        C_matrix_1x2_int_type,
        B_matrix_2x1_int_type
    >::type CxB;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<int>
        >
    > CxB_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<CxB, CxB_manual_product>));

    typedef bub::result_of::matrix_product<
        B_matrix_2x1_int_type,
        D_matrix_1x2_int_double_type
    >::type BxD;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<int, double>,
            boost::fusion::vector<int, double>
        >
    > BxD_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<BxD, BxD_manual_product>));

    typedef bub::result_of::matrix_product<
        D_matrix_1x2_int_double_type,
        B_matrix_2x1_int_type
    >::type DxB;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<double>
        >
    > DxB_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<DxB, DxB_manual_product>));

    typedef bub::result_of::matrix_product<
        C_matrix_1x2_int_type,
        E_matrix_2x2_double_type
    >::type CxE;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<double, double>
        >
    > CxE_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<CxE, CxE_manual_product>));

    typedef bub::result_of::matrix_product<
        E_matrix_2x2_double_type,
        B_matrix_2x1_int_type
    >::type ExB;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<double>,
            boost::fusion::vector<double>
        >
    > ExB_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<ExB, ExB_manual_product>));

    typedef bub::result_of::matrix_product<
        E_matrix_2x2_double_type,
        E_matrix_2x2_double_type
    >::type ExE;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<double, double>,
            boost::fusion::vector<double, double>
        >
    > ExE_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<ExE, ExE_manual_product>));


    // matrix-scalar products

    typedef bub::result_of::scalar_product<
        A_matrix_2x1_float_int_type,
        int
    >::type Axint;
    BOOST_MPL_ASSERT((boost::is_same<Axint, A_matrix_2x1_float_int_type>));

    typedef bub::result_of::scalar_product<
        A_matrix_2x1_float_int_type,
        double
    >::type Axdouble;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<double>,
            boost::fusion::vector<double>
        >
    > Axdouble_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<Axdouble, Axdouble_manual_product>));

    typedef bub::result_of::scalar_product<
        D_matrix_1x2_int_double_type,
        int
    >::type Dxint;
    BOOST_MPL_ASSERT((boost::is_same<Dxint, D_matrix_1x2_int_double_type>));

    typedef bub::result_of::scalar_product<
        D_matrix_1x2_int_double_type,
        double
    >::type Dxdouble;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<double, double>
        >
    > Dxdouble_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<Dxdouble, Dxdouble_manual_product>));


    // matrix-scalar quotients

    typedef bub::result_of::scalar_quotient<
        A_matrix_2x1_float_int_type,
        int
    >::type Adint;
    BOOST_MPL_ASSERT((boost::is_same<Adint, A_matrix_2x1_float_int_type>));

    typedef bub::result_of::scalar_quotient<
        A_matrix_2x1_float_int_type,
        double
    >::type Addouble;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<double>,
            boost::fusion::vector<double>
        >
    > Addouble_manual_quotient;
    BOOST_MPL_ASSERT((boost::is_same<Addouble, Addouble_manual_quotient>));

    typedef bub::result_of::scalar_quotient<
        D_matrix_1x2_int_double_type,
        int
    >::type Ddint;
    BOOST_MPL_ASSERT((boost::is_same<Ddint, D_matrix_1x2_int_double_type>));

    typedef bub::result_of::scalar_quotient<
        D_matrix_1x2_int_double_type,
        double
    >::type Dddouble;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<double, double>
        >
    > Dddouble_manual_quotient;
    BOOST_MPL_ASSERT((boost::is_same<Dddouble, Dddouble_manual_quotient>));


    // element-wise matrix sums

    typedef bub::result_of::matrix_element_sum<
        A_matrix_2x1_float_int_type,
        B_matrix_2x1_int_type
    >::type ApB;
    BOOST_MPL_ASSERT((boost::is_same<ApB, A_matrix_2x1_float_int_type>));

    typedef bub::result_of::matrix_element_sum<
        C_matrix_1x2_int_type,
        D_matrix_1x2_int_double_type
    >::type CpD;
    BOOST_MPL_ASSERT((boost::is_same<CpD, D_matrix_1x2_int_double_type>));


    // element-wise matrix differences

    typedef bub::result_of::matrix_element_difference<
        A_matrix_2x1_float_int_type,
        B_matrix_2x1_int_type
    >::type AmB;
    BOOST_MPL_ASSERT((boost::is_same<ApB, A_matrix_2x1_float_int_type>));

    typedef bub::result_of::matrix_element_difference<
        C_matrix_1x2_int_type,
        D_matrix_1x2_int_double_type
    >::type CmD;
    BOOST_MPL_ASSERT((boost::is_same<CpD, D_matrix_1x2_int_double_type>));


    // element-wise matrix products

    typedef bub::result_of::matrix_element_product<
        A_matrix_2x1_float_int_type,
        B_matrix_2x1_int_type
    >::type AxB_elementwise_product;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<float>,
            boost::fusion::vector<int>
        >
    > AxB_elementwise_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<AxB_elementwise_product, AxB_elementwise_manual_product>));

    typedef bub::result_of::matrix_element_product<
        C_matrix_1x2_int_type,
        D_matrix_1x2_int_double_type
    >::type CxD_elementwise_product;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<int, double>
        >
    > CxD_elementwise_manual_product;
    BOOST_MPL_ASSERT((boost::is_same<CxD_elementwise_product, CxD_elementwise_manual_product>));

    // element-wise matrix quotients

    typedef bub::result_of::matrix_element_quotient<
        A_matrix_2x1_float_int_type,
        B_matrix_2x1_int_type
    >::type AxB_elementwise_quotient;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<float>,
            boost::fusion::vector<int>
        >
    > AxB_elementwise_manual_quotient;
    BOOST_MPL_ASSERT((boost::is_same<AxB_elementwise_quotient, AxB_elementwise_manual_quotient>));

    typedef bub::result_of::matrix_element_quotient<
        C_matrix_1x2_int_type,
        D_matrix_1x2_int_double_type
    >::type CxD_elementwise_quotient;
    typedef bub::matrix<
        boost::fusion::vector<
            boost::fusion::vector<int, double>
        >
    > CxD_elementwise_manual_quotient;
    BOOST_MPL_ASSERT((boost::is_same<CxD_elementwise_quotient, CxD_elementwise_manual_quotient>));

    return 0;
}
