// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units/is_dimensionless.hpp>
#include <boost/units_blas/identity_matrix.hpp>
#include <boost/units_blas/operations.hpp>
#include <boost/units_blas/traits.hpp>
#include <boost/units_blas/result_of/at.hpp>

#include <boost/test/minimal.hpp>


typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, long, int>,
        boost::fusion::vector<long, float, int>,
        boost::fusion::vector<float, float, long>
    >
>::type matrix_3x3_type;

typedef boost::units::quantity<boost::units::si::length> length;
typedef boost::units::quantity<boost::units::si::velocity> velocity;
typedef boost::units::quantity<boost::units::si::dimensionless> dimensionless;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, length>,
        boost::fusion::vector<length, length>
    >
>::type matrix_2x2_units_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<length_sq, length_sq, bub::_               >,
        boost::fusion::vector<length_sq, length_sq, bub::_               >,
        boost::fusion::vector<bub::_,    bub::_,    length_sq_per_time_sq>
    >
>::type matrix_3x3_units_type;

int test_main (int, char *[])
{
    {
        matrix_3x3_type matrix;
        matrix.at<0, 0>() = 1;
        matrix.at<0, 1>() = 2;
        matrix.at<0, 2>() = 3;
        matrix.at<1, 0>() = 4;
        matrix.at<1, 1>() = 5.0f;
        matrix.at<1, 2>() = 6;
        matrix.at<2, 0>() = 7.0f;
        matrix.at<2, 1>() = 8.0f;
        matrix.at<2, 2>() = 9;

        typedef BOOST_TYPEOF((bub::make_identity_matrix<matrix_3x3_type>())) identity_type;
        identity_type identity_matrix = bub::make_identity_matrix<matrix_3x3_type>();

        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<0, 0>())), long>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<0, 1>())), float>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<0, 2>())), float>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<1, 0>())), float>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<1, 1>())), float>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<1, 2>())), float>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<2, 0>())), float>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<2, 1>())), float>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<2, 2>())), float>));

        BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
        BOOST_CHECK((identity_matrix.at<0, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<0, 2>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 1>() == 1));
        BOOST_CHECK((identity_matrix.at<1, 2>() == 0));
        BOOST_CHECK((identity_matrix.at<2, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<2, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<2, 2>() == 1));

        matrix_3x3_type matrix_times_identity = matrix * identity_matrix;

        BOOST_CHECK((matrix_times_identity.at<0, 0>() == 1));
        BOOST_CHECK((matrix_times_identity.at<0, 1>() == 2));
        BOOST_CHECK((matrix_times_identity.at<0, 2>() == 3));
        BOOST_CHECK((matrix_times_identity.at<1, 0>() == 4));
        BOOST_CHECK((matrix_times_identity.at<1, 1>() == 5));
        BOOST_CHECK((matrix_times_identity.at<1, 2>() == 6));
        BOOST_CHECK((matrix_times_identity.at<2, 0>() == 7));
        BOOST_CHECK((matrix_times_identity.at<2, 1>() == 8));
        BOOST_CHECK((matrix_times_identity.at<2, 2>() == 9));
    }

    {
        matrix_2x2_units_type matrix;
        matrix.at<0, 0>() = length::from_value(1.0);
        matrix.at<0, 1>() = length::from_value(2.0);
        matrix.at<1, 0>() = length::from_value(3.0);
        matrix.at<1, 1>() = length::from_value(4.0);

        typedef BOOST_TYPEOF((bub::make_identity_matrix<matrix_2x2_units_type>())) identity_type;
        identity_type identity_matrix = bub::make_identity_matrix<matrix_2x2_units_type>();

        BOOST_MPL_ASSERT((boost::units::is_dimensionless<BOOST_TYPEOF((identity_matrix.at<0, 0>()))>));
        BOOST_MPL_ASSERT((boost::units::is_dimensionless<BOOST_TYPEOF((identity_matrix.at<0, 1>()))>));
        BOOST_MPL_ASSERT((boost::units::is_dimensionless<BOOST_TYPEOF((identity_matrix.at<1, 0>()))>));
        BOOST_MPL_ASSERT((boost::units::is_dimensionless<BOOST_TYPEOF((identity_matrix.at<1, 1>()))>));

        BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
        BOOST_CHECK((identity_matrix.at<0, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 1>() == 1));

        matrix_2x2_units_type matrix_times_identity = matrix * identity_matrix;

        BOOST_CHECK((matrix_times_identity.at<0, 0>() == length::from_value(1.0)));
        BOOST_CHECK((matrix_times_identity.at<0, 1>() == length::from_value(2.0)));
        BOOST_CHECK((matrix_times_identity.at<1, 0>() == length::from_value(3.0)));
        BOOST_CHECK((matrix_times_identity.at<1, 1>() == length::from_value(4.0)));
    }

    {
        matrix_3x3_units_type matrix;
        matrix.at<0, 0>() = length_sq::from_value(1.0);
        matrix.at<0, 1>() = length_sq::from_value(2.0);
        matrix.at<1, 0>() = length_sq::from_value(3.0);
        matrix.at<1, 1>() = length_sq::from_value(4.0);
        matrix.at<2, 2>() = length_sq_per_time_sq::from_value(5.0);

        typedef BOOST_TYPEOF((bub::make_identity_matrix<matrix_3x3_units_type>())) identity_type;
        identity_type identity_matrix = bub::make_identity_matrix<matrix_3x3_units_type>();

        BOOST_MPL_ASSERT((boost::units::is_dimensionless<BOOST_TYPEOF((identity_matrix.at<0, 0>()))>));
        BOOST_MPL_ASSERT((boost::units::is_dimensionless<BOOST_TYPEOF((identity_matrix.at<0, 1>()))>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<0, 2>())), bub::_>));
        BOOST_MPL_ASSERT((boost::units::is_dimensionless<BOOST_TYPEOF((identity_matrix.at<1, 0>()))>));
        BOOST_MPL_ASSERT((boost::units::is_dimensionless<BOOST_TYPEOF((identity_matrix.at<1, 1>()))>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<1, 2>())), bub::_>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<2, 0>())), bub::_>));
        BOOST_MPL_ASSERT((boost::is_same<BOOST_TYPEOF((identity_matrix.at<2, 1>())), bub::_>));
        BOOST_MPL_ASSERT((boost::units::is_dimensionless<BOOST_TYPEOF((identity_matrix.at<2, 2>()))>));

        BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
        BOOST_CHECK((identity_matrix.at<0, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 1>() == 1));
        BOOST_CHECK((identity_matrix.at<2, 2>() == 1));

        matrix_3x3_units_type matrix_times_identity = matrix * identity_matrix;

        BOOST_CHECK((matrix_times_identity.at<0, 0>() == length_sq::from_value(1.0)));
        BOOST_CHECK((matrix_times_identity.at<0, 1>() == length_sq::from_value(2.0)));
        BOOST_CHECK((matrix_times_identity.at<1, 0>() == length_sq::from_value(3.0)));
        BOOST_CHECK((matrix_times_identity.at<1, 1>() == length_sq::from_value(4.0)));
        BOOST_CHECK((matrix_times_identity.at<2, 2>() == length_sq_per_time_sq::from_value(5.0)));
    }

    return 0;
}
