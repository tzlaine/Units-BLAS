// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/traits.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

#include <boost/test/minimal.hpp>


namespace bub = boost::units_blas;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, long, float>,
        boost::fusion::vector<long, float, int>,
        boost::fusion::vector<int, float, long>
    >
> matrix_3x3_type;

typedef boost::units::quantity<boost::units::SI::length> length;
typedef boost::units::quantity<boost::units::SI::velocity> velocity;
typedef boost::units::quantity<boost::units::SI::dimensionless> dimensionless;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, velocity, dimensionless>,
        boost::fusion::vector<velocity, dimensionless, length>,
        boost::fusion::vector<length, dimensionless, velocity>
    >
> matrix_3x3_units_type;

int test_main (int, char *[])
{
    matrix_3x3_type identity_matrix = bub::identity_matrix<matrix_3x3_type>();

    typedef BOOST_TYPEOF((bub::identity_matrix<matrix_3x3_type>())) identity_type;
    BOOST_MPL_ASSERT((boost::is_same<matrix_3x3_type, identity_type>));

    BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
    BOOST_CHECK((identity_matrix.at<0, 1>() == 0l));
    BOOST_CHECK((identity_matrix.at<0, 2>() == 0.0f));
    BOOST_CHECK((identity_matrix.at<1, 0>() == 0l));
    BOOST_CHECK((identity_matrix.at<1, 1>() == 1.0));
    BOOST_CHECK((identity_matrix.at<1, 2>() == 0));
    BOOST_CHECK((identity_matrix.at<2, 0>() == 0));
    BOOST_CHECK((identity_matrix.at<2, 1>() == 0.0f));
    BOOST_CHECK((identity_matrix.at<2, 2>() == 1l));

    matrix_3x3_units_type identity_units_matrix = bub::identity_matrix<matrix_3x3_units_type>();

    typedef BOOST_TYPEOF((bub::identity_matrix<matrix_3x3_units_type>())) identity_units_type;
    BOOST_MPL_ASSERT((boost::is_same<matrix_3x3_units_type, identity_units_type>));

    BOOST_CHECK((identity_units_matrix.at<0, 0>() == length::from_value(1.0)));
    BOOST_CHECK((identity_units_matrix.at<0, 1>() == velocity::from_value(0.0)));
    BOOST_CHECK((identity_units_matrix.at<0, 2>() == dimensionless::from_value(0.0)));
    BOOST_CHECK((identity_units_matrix.at<1, 0>() == velocity::from_value(0.0)));
    BOOST_CHECK((identity_units_matrix.at<1, 1>() == dimensionless::from_value(1.0)));
    BOOST_CHECK((identity_units_matrix.at<1, 2>() == length::from_value(0.0)));
    BOOST_CHECK((identity_units_matrix.at<2, 0>() == length::from_value(0.0)));
    BOOST_CHECK((identity_units_matrix.at<2, 1>() == dimensionless::from_value(0.0)));
    BOOST_CHECK((identity_units_matrix.at<2, 2>() == velocity::from_value(1.0)));

    return 0;
}
