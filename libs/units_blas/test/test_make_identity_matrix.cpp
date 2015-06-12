// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units/is_dimensionless.hpp>
#include <boost/units_blas/identity_matrix.hpp>
#include <boost/units_blas/operations.hpp>

#include <boost/test/minimal.hpp>


using matrix_3x3_type = bub::matrix<
    bh::_tuple<int, long, int>,
    bh::_tuple<long, float, int>,
    bh::_tuple<float, float, long>
>;

using matrix_2x2_units_type = bub::matrix<
    bh::_tuple<length, length>,
    bh::_tuple<length, length>
>;

using xyd_matrix_3x3_units_type = bub::matrix<
    bh::_tuple<length_sq,          length_sq,          length_sq_per_time>,
    bh::_tuple<length_sq,          length_sq,          length_sq_per_time>,
    bh::_tuple<length_sq_per_time, length_sq_per_time, length_sq_per_time_sq>
>;

#define CHECK_TYPE(i, j, T)                                             \
        do {                                                            \
            auto x = identity_matrix.at<i, j>();                        \
            static_assert(std::is_same<decltype(x), T>::value, "");     \
        } while (false)

#define CHECK_DIMENSIONLESS(i, j)                                       \
        do {                                                            \
            auto x = identity_matrix.at<i, j>();                        \
            static_assert(boost::units::is_dimensionless<decltype(x)>::value, ""); \
        } while (false)

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

        auto identity_matrix = bub::make_left_identity_matrix<matrix_3x3_type>();

        CHECK_TYPE(0, 0, long);
        CHECK_TYPE(0, 1, float);
        CHECK_TYPE(0, 2, float);
        CHECK_TYPE(1, 0, float);
        CHECK_TYPE(1, 1, float);
        CHECK_TYPE(1, 2, float);
        CHECK_TYPE(2, 0, float);
        CHECK_TYPE(2, 1, float);
        CHECK_TYPE(2, 2, float);

        BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
        BOOST_CHECK((identity_matrix.at<0, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<0, 2>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 1>() == 1));
        BOOST_CHECK((identity_matrix.at<1, 2>() == 0));
        BOOST_CHECK((identity_matrix.at<2, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<2, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<2, 2>() == 1));

        matrix_3x3_type matrix_times_left_identity = matrix * identity_matrix;

        BOOST_CHECK((matrix_times_left_identity.at<0, 0>() == 1));
        BOOST_CHECK((matrix_times_left_identity.at<0, 1>() == 2));
        BOOST_CHECK((matrix_times_left_identity.at<0, 2>() == 3));
        BOOST_CHECK((matrix_times_left_identity.at<1, 0>() == 4));
        BOOST_CHECK((matrix_times_left_identity.at<1, 1>() == 5));
        BOOST_CHECK((matrix_times_left_identity.at<1, 2>() == 6));
        BOOST_CHECK((matrix_times_left_identity.at<2, 0>() == 7));
        BOOST_CHECK((matrix_times_left_identity.at<2, 1>() == 8));
        BOOST_CHECK((matrix_times_left_identity.at<2, 2>() == 9));
    }

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

        auto identity_matrix = bub::make_right_identity_matrix<matrix_3x3_type>();

        CHECK_TYPE(0, 0, float);
        CHECK_TYPE(0, 1, float);
        CHECK_TYPE(0, 2, float);
        CHECK_TYPE(1, 0, float);
        CHECK_TYPE(1, 1, float);
        CHECK_TYPE(1, 2, float);
        CHECK_TYPE(2, 0, float);
        CHECK_TYPE(2, 1, float);
        CHECK_TYPE(2, 2, long);

        BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
        BOOST_CHECK((identity_matrix.at<0, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<0, 2>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 1>() == 1));
        BOOST_CHECK((identity_matrix.at<1, 2>() == 0));
        BOOST_CHECK((identity_matrix.at<2, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<2, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<2, 2>() == 1));

        matrix_3x3_type matrix_times_right_identity = matrix * identity_matrix;

        BOOST_CHECK((matrix_times_right_identity.at<0, 0>() == 1));
        BOOST_CHECK((matrix_times_right_identity.at<0, 1>() == 2));
        BOOST_CHECK((matrix_times_right_identity.at<0, 2>() == 3));
        BOOST_CHECK((matrix_times_right_identity.at<1, 0>() == 4));
        BOOST_CHECK((matrix_times_right_identity.at<1, 1>() == 5));
        BOOST_CHECK((matrix_times_right_identity.at<1, 2>() == 6));
        BOOST_CHECK((matrix_times_right_identity.at<2, 0>() == 7));
        BOOST_CHECK((matrix_times_right_identity.at<2, 1>() == 8));
        BOOST_CHECK((matrix_times_right_identity.at<2, 2>() == 9));
    }

    {
        matrix_2x2_units_type matrix;
        matrix.at<0, 0>() = length::from_value(1.0);
        matrix.at<0, 1>() = length::from_value(2.0);
        matrix.at<1, 0>() = length::from_value(3.0);
        matrix.at<1, 1>() = length::from_value(4.0);

        auto identity_matrix = bub::make_left_identity_matrix<matrix_2x2_units_type>();

        CHECK_DIMENSIONLESS(0, 0);
        CHECK_DIMENSIONLESS(0, 1);
        CHECK_DIMENSIONLESS(1, 0);
        CHECK_DIMENSIONLESS(1, 1);

        BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
        BOOST_CHECK((identity_matrix.at<0, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 1>() == 1));

        matrix_2x2_units_type matrix_times_left_identity = matrix * identity_matrix;

        BOOST_CHECK((matrix_times_left_identity.at<0, 0>() == length::from_value(1.0)));
        BOOST_CHECK((matrix_times_left_identity.at<0, 1>() == length::from_value(2.0)));
        BOOST_CHECK((matrix_times_left_identity.at<1, 0>() == length::from_value(3.0)));
        BOOST_CHECK((matrix_times_left_identity.at<1, 1>() == length::from_value(4.0)));
    }

    {
        matrix_2x2_units_type matrix;
        matrix.at<0, 0>() = length::from_value(1.0);
        matrix.at<0, 1>() = length::from_value(2.0);
        matrix.at<1, 0>() = length::from_value(3.0);
        matrix.at<1, 1>() = length::from_value(4.0);

        auto identity_matrix = bub::make_right_identity_matrix<matrix_2x2_units_type>();

        CHECK_DIMENSIONLESS(0, 0);
        CHECK_DIMENSIONLESS(0, 1);
        CHECK_DIMENSIONLESS(1, 0);
        CHECK_DIMENSIONLESS(1, 1);

        BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
        BOOST_CHECK((identity_matrix.at<0, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 1>() == 1));

        matrix_2x2_units_type matrix_times_right_identity = matrix * identity_matrix;

        BOOST_CHECK((matrix_times_right_identity.at<0, 0>() == length::from_value(1.0)));
        BOOST_CHECK((matrix_times_right_identity.at<0, 1>() == length::from_value(2.0)));
        BOOST_CHECK((matrix_times_right_identity.at<1, 0>() == length::from_value(3.0)));
        BOOST_CHECK((matrix_times_right_identity.at<1, 1>() == length::from_value(4.0)));
    }

    {
        xyd_matrix_3x3_units_type matrix;
        matrix.at<0, 0>() = length_sq::from_value(1.0);
        matrix.at<0, 1>() = length_sq::from_value(2.0);
        matrix.at<0, 2>() = length_sq_per_time::from_value(3.0);
        matrix.at<1, 0>() = length_sq::from_value(4.0);
        matrix.at<1, 1>() = length_sq::from_value(5.0);
        matrix.at<1, 2>() = length_sq_per_time::from_value(6.0);
        matrix.at<2, 0>() = length_sq_per_time::from_value(7.0);
        matrix.at<2, 1>() = length_sq_per_time::from_value(8.0);
        matrix.at<2, 2>() = length_sq_per_time_sq::from_value(9.0);

        auto identity_matrix = bub::make_left_identity_matrix<xyd_matrix_3x3_units_type>();

        CHECK_DIMENSIONLESS(0, 0);
        CHECK_DIMENSIONLESS(0, 1);
        CHECK_TYPE(0, 2, time_);
        CHECK_DIMENSIONLESS(1, 0);
        CHECK_DIMENSIONLESS(1, 1);
        CHECK_TYPE(1, 2, time_);
        CHECK_TYPE(2, 0, frequency);
        CHECK_TYPE(2, 1, frequency);
        CHECK_DIMENSIONLESS(2, 2);

        BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
        BOOST_CHECK((identity_matrix.at<0, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<0, 2>() == time_::from_value(0)));
        BOOST_CHECK((identity_matrix.at<1, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 1>() == 1));
        BOOST_CHECK((identity_matrix.at<1, 2>() == time_::from_value(0)));
        BOOST_CHECK((identity_matrix.at<2, 0>() == frequency::from_value(0)));
        BOOST_CHECK((identity_matrix.at<2, 1>() == frequency::from_value(0)));
        BOOST_CHECK((identity_matrix.at<2, 2>() == 1));

        xyd_matrix_3x3_units_type matrix_times_left_identity = identity_matrix * matrix;

        BOOST_CHECK((matrix_times_left_identity.at<0, 0>() == length_sq::from_value(1.0)));
        BOOST_CHECK((matrix_times_left_identity.at<0, 1>() == length_sq::from_value(2.0)));
        BOOST_CHECK((matrix_times_left_identity.at<0, 2>() == length_sq_per_time::from_value(3.0)));
        BOOST_CHECK((matrix_times_left_identity.at<1, 0>() == length_sq::from_value(4.0)));
        BOOST_CHECK((matrix_times_left_identity.at<1, 1>() == length_sq::from_value(5.0)));
        BOOST_CHECK((matrix_times_left_identity.at<1, 2>() == length_sq_per_time::from_value(6.0)));
        BOOST_CHECK((matrix_times_left_identity.at<2, 0>() == length_sq_per_time::from_value(7.0)));
        BOOST_CHECK((matrix_times_left_identity.at<2, 1>() == length_sq_per_time::from_value(8.0)));
        BOOST_CHECK((matrix_times_left_identity.at<2, 2>() == length_sq_per_time_sq::from_value(9.0)));
    }

    {
        xyd_matrix_3x3_units_type matrix;
        matrix.at<0, 0>() = length_sq::from_value(1.0);
        matrix.at<0, 1>() = length_sq::from_value(2.0);
        matrix.at<0, 2>() = length_sq_per_time::from_value(3.0);
        matrix.at<1, 0>() = length_sq::from_value(4.0);
        matrix.at<1, 1>() = length_sq::from_value(5.0);
        matrix.at<1, 2>() = length_sq_per_time::from_value(6.0);
        matrix.at<2, 0>() = length_sq_per_time::from_value(7.0);
        matrix.at<2, 1>() = length_sq_per_time::from_value(8.0);
        matrix.at<2, 2>() = length_sq_per_time_sq::from_value(9.0);

        auto identity_matrix = bub::make_right_identity_matrix<xyd_matrix_3x3_units_type>();

        CHECK_DIMENSIONLESS(0, 0);
        CHECK_DIMENSIONLESS(0, 1);
        CHECK_TYPE(0, 2, frequency);
        CHECK_DIMENSIONLESS(1, 0);
        CHECK_DIMENSIONLESS(1, 1);
        CHECK_TYPE(1, 2, frequency);
        CHECK_TYPE(2, 0, time_);
        CHECK_TYPE(2, 1, time_);
        CHECK_DIMENSIONLESS(2, 2);

        BOOST_CHECK((identity_matrix.at<0, 0>() == 1));
        BOOST_CHECK((identity_matrix.at<0, 1>() == 0));
        BOOST_CHECK((identity_matrix.at<0, 2>() == frequency::from_value(0)));
        BOOST_CHECK((identity_matrix.at<1, 0>() == 0));
        BOOST_CHECK((identity_matrix.at<1, 1>() == 1));
        BOOST_CHECK((identity_matrix.at<1, 2>() == frequency::from_value(0)));
        BOOST_CHECK((identity_matrix.at<2, 0>() == time_::from_value(0)));
        BOOST_CHECK((identity_matrix.at<2, 1>() == time_::from_value(0)));
        BOOST_CHECK((identity_matrix.at<2, 2>() == 1));

        xyd_matrix_3x3_units_type matrix_times_right_identity = matrix * identity_matrix;

        BOOST_CHECK((matrix_times_right_identity.at<0, 0>() == length_sq::from_value(1.0)));
        BOOST_CHECK((matrix_times_right_identity.at<0, 1>() == length_sq::from_value(2.0)));
        BOOST_CHECK((matrix_times_right_identity.at<0, 2>() == length_sq_per_time::from_value(3.0)));
        BOOST_CHECK((matrix_times_right_identity.at<1, 0>() == length_sq::from_value(4.0)));
        BOOST_CHECK((matrix_times_right_identity.at<1, 1>() == length_sq::from_value(5.0)));
        BOOST_CHECK((matrix_times_right_identity.at<1, 2>() == length_sq_per_time::from_value(6.0)));
        BOOST_CHECK((matrix_times_right_identity.at<2, 0>() == length_sq_per_time::from_value(7.0)));
        BOOST_CHECK((matrix_times_right_identity.at<2, 1>() == length_sq_per_time::from_value(8.0)));
        BOOST_CHECK((matrix_times_right_identity.at<2, 2>() == length_sq_per_time_sq::from_value(9.0)));
    }

    return 0;
}
