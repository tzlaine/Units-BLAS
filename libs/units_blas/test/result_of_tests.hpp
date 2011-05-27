// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_UNITS_BLAS_TEST_RESULT_OF_TESTS_HPP
#define BOOST_UNITS_BLAS_TEST_RESULT_OF_TESTS_HPP

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/make_matrix.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/frequency.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

namespace bub = boost::units_blas;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<float>,
        boost::fusion::vector<int>
    >
>::type A_matrix_2x1_float_int_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<int>,
        boost::fusion::vector<int>
    >
>::type B_matrix_2x1_int_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, int>
    >
>::type C_matrix_1x2_int_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, double>
    >
>::type D_matrix_1x2_int_double_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<double, double>,
        boost::fusion::vector<double, double>
    >
>::type E_matrix_2x2_double_type;

typedef boost::units::quantity<boost::units::si::time> time_;
typedef boost::units::quantity<boost::units::si::length> length;
typedef boost::units::quantity<boost::units::si::velocity> velocity;
typedef boost::units::quantity<boost::units::si::dimensionless> dimensionless;
typedef boost::units::quantity<boost::units::si::frequency> frequency;

namespace boost { namespace units { namespace si {
    typedef derived_dimension<time_base_dimension, 1, length_base_dimension, 1>::type time_length_dimension;
    typedef derived_dimension<time_base_dimension, 2>::type time_sq_dimension;
    typedef derived_dimension<length_base_dimension, 2>::type length_sq_dimension;
    typedef derived_dimension<time_base_dimension, -1, length_base_dimension, 2>::type length_sq_per_time_dimension;
    typedef derived_dimension<time_base_dimension, -2, length_base_dimension, 2>::type length_sq_per_time_sq_dimension;
    typedef derived_dimension<length_base_dimension, -2>::type inv_length_sq_dimension;
    typedef derived_dimension<time_base_dimension, 1, length_base_dimension, -2>::type time_per_length_sq_dimension;
    typedef derived_dimension<time_base_dimension, 2, length_base_dimension, -2>::type time_sq_per_length_sq_dimension;

    typedef unit<time_length_dimension, system> time_length;
    typedef unit<time_sq_dimension, system> time_sq;
    typedef unit<length_sq_dimension, system> length_sq;
    typedef unit<length_sq_per_time_dimension, system> length_sq_per_time;
    typedef unit<length_sq_per_time_sq_dimension, system> length_sq_per_time_sq;
    typedef unit<inv_length_sq_dimension, system> inv_length_sq;
    typedef unit<time_per_length_sq_dimension, system> time_per_length_sq;
    typedef unit<time_sq_per_length_sq_dimension, system> time_sq_per_length_sq;
} } }

typedef boost::units::quantity<boost::units::si::time_length> time_length;
typedef boost::units::quantity<boost::units::si::time_sq> time_sq;
typedef boost::units::quantity<boost::units::si::length_sq> length_sq;
typedef boost::units::quantity<boost::units::si::length_sq_per_time> length_sq_per_time;
typedef boost::units::quantity<boost::units::si::length_sq_per_time_sq> length_sq_per_time_sq;
typedef boost::units::quantity<boost::units::si::inv_length_sq> inv_length_sq;
typedef boost::units::quantity<boost::units::si::time_per_length_sq> time_per_length_sq;
typedef boost::units::quantity<boost::units::si::time_sq_per_length_sq> time_sq_per_length_sq;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_>,
        boost::fusion::vector<length>
    >
>::type A_matrix_2x1_time_length_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<length>,
        boost::fusion::vector<length>
    >
>::type B_matrix_2x1_length_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, length>
    >
>::type C_matrix_1x2_length_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, time_>
    >
>::type D_matrix_1x2_length_time_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_, time_>,
        boost::fusion::vector<time_, time_>
    >
>::type E_matrix_2x2_time_type;

#endif // BOOST_UNITS_BLAS_TEST_RESULT_OF_TESTS_HPP
