// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_UNITS_BLAS_TEST_OPERATIONS_TESTS_HPP
#define BOOST_UNITS_BLAS_TEST_OPERATIONS_TESTS_HPP

#define BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS 1

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/make_matrix.hpp>
#include <boost/units_blas/operations.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/frequency.hpp>
#include <boost/units/systems/si/dimensionless.hpp>


namespace bub = boost::units_blas;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<double>,
        boost::fusion::vector<double>,
        boost::fusion::vector<double>
    >
>::type A_matrix_3x1_fundamentals_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<double, double, double>
    >
>::type B_matrix_1x3_fundamentals_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<double, double, double, double>,
        boost::fusion::vector<double, double, double, double>,
        boost::fusion::vector<double, double, double, double>,
        boost::fusion::vector<double, double, double, double>
    >
>::type D_matrix_4x4_fundamentals_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<double>
    >
>::type E_matrix_1x1_fundamentals_type;

typedef boost::units::quantity<boost::units::si::time> time_;
typedef boost::units::quantity<boost::units::si::length> length;
typedef boost::units::quantity<boost::units::si::velocity> velocity;
typedef boost::units::quantity<boost::units::si::dimensionless> dimensionless;
typedef boost::units::quantity<boost::units::si::frequency> frequency;

namespace boost { namespace units { namespace si {
    typedef derived_dimension<time_base_dimension, 1, length_base_dimension, 1>::type time_length_dimension;
    typedef derived_dimension<time_base_dimension, 2>::type time_sq_dimension;
    typedef derived_dimension<length_base_dimension, 2>::type length_sq_dimension;
    typedef unit<time_length_dimension, system> time_length;
    typedef unit<time_sq_dimension, system> time_sq;
    typedef unit<length_sq_dimension, system> length_sq;
} } }

typedef boost::units::quantity<boost::units::si::time_length> time_length;
typedef boost::units::quantity<boost::units::si::time_sq> time_sq;
typedef boost::units::quantity<boost::units::si::length_sq> length_sq;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_>,
        boost::fusion::vector<length>,
        boost::fusion::vector<dimensionless>
    >
>::type A_matrix_3x1_units_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<length>,
        boost::fusion::vector<time_>,
        boost::fusion::vector<time_length>
    >
>::type A_matrix_3x1_units_type_2;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_, time_, time_>
    >
>::type B_matrix_1x3_units_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<dimensionless>,
        boost::fusion::vector<dimensionless>,
        boost::fusion::vector<dimensionless>
    >
>::type C_matrix_3x1_units_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, length, length, length>,
        boost::fusion::vector<length, length, length, length>,
        boost::fusion::vector<length, length, length, length>,
        boost::fusion::vector<length, length, length, length>
    >
>::type D_matrix_4x4_units_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<length>
    >
>::type E_matrix_1x1_units_type;

#endif // BOOST_UNITS_BLAS_TEST_OPERATIONS_TESTS_HPP
