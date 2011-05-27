// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "operations_tests.hpp"

#include <boost/test/minimal.hpp>


int test_main (int, char *[])
{
    // fundamental types

    A_matrix_3x1_fundamentals_type m_A1;
    m_A1.at<0, 0>() = 1.0;
    m_A1.at<1, 0>() = 2.0;
    m_A1.at<2, 0>() = 3.0;
    A_matrix_3x1_fundamentals_type m_A2 = m_A1;
    BOOST_CHECK((dot(m_A1, m_A2) == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));
    BOOST_CHECK((m_A1 * m_A2 == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));

    B_matrix_1x3_fundamentals_type m_B1;
    m_B1.at<0, 0>() = 1.0;
    m_B1.at<0, 1>() = 2.0;
    m_B1.at<0, 2>() = 3.0;
    B_matrix_1x3_fundamentals_type m_B2 = m_B1;
    BOOST_CHECK((dot(m_B1, m_B2) == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));
    BOOST_CHECK((m_B1 * m_B2 == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));

    A_matrix_3x1_fundamentals_type m_A_x_dir;
    m_A_x_dir.at<0, 0>() = 1.0;
    m_A_x_dir.at<1, 0>() = 0.0;
    m_A_x_dir.at<2, 0>() = 0.0;
    A_matrix_3x1_fundamentals_type m_A_y_dir;
    m_A_y_dir.at<0, 0>() = 0.0;
    m_A_y_dir.at<1, 0>() = 1.0;
    m_A_y_dir.at<2, 0>() = 0.0;
    A_matrix_3x1_fundamentals_type m_A_z_dir;
    m_A_z_dir.at<0, 0>() = 0.0;
    m_A_z_dir.at<1, 0>() = 0.0;
    m_A_z_dir.at<2, 0>() = 1.0;

    A_matrix_3x1_fundamentals_type m_A_x_cross_y = cross(m_A_x_dir, m_A_y_dir);
    BOOST_CHECK((m_A_x_cross_y.at<0, 0>() == 0.0));
    BOOST_CHECK((m_A_x_cross_y.at<1, 0>() == 0.0));
    BOOST_CHECK((m_A_x_cross_y.at<2, 0>() == 1.0));
    m_A_x_cross_y = m_A_x_dir ^ m_A_y_dir;
    BOOST_CHECK((m_A_x_cross_y.at<0, 0>() == 0.0));
    BOOST_CHECK((m_A_x_cross_y.at<1, 0>() == 0.0));
    BOOST_CHECK((m_A_x_cross_y.at<2, 0>() == 1.0));
    A_matrix_3x1_fundamentals_type m_A_y_cross_x = cross(m_A_y_dir, m_A_x_dir);
    BOOST_CHECK((m_A_y_cross_x.at<0, 0>() == 0.0));
    BOOST_CHECK((m_A_y_cross_x.at<1, 0>() == 0.0));
    BOOST_CHECK((m_A_y_cross_x.at<2, 0>() == -1.0));
    m_A_y_cross_x = m_A_y_dir ^ m_A_x_dir;
    BOOST_CHECK((m_A_y_cross_x.at<0, 0>() == 0.0));
    BOOST_CHECK((m_A_y_cross_x.at<1, 0>() == 0.0));
    BOOST_CHECK((m_A_y_cross_x.at<2, 0>() == -1.0));
    A_matrix_3x1_fundamentals_type m_A_x_cross_z = cross(m_A_x_dir, m_A_z_dir);
    BOOST_CHECK((m_A_x_cross_z.at<0, 0>() == 0.0));
    BOOST_CHECK((m_A_x_cross_z.at<1, 0>() == -1.0));
    BOOST_CHECK((m_A_x_cross_z.at<2, 0>() == 0.0));
    m_A_x_cross_z = m_A_x_dir ^ m_A_z_dir;
    BOOST_CHECK((m_A_x_cross_z.at<0, 0>() == 0.0));
    BOOST_CHECK((m_A_x_cross_z.at<1, 0>() == -1.0));
    BOOST_CHECK((m_A_x_cross_z.at<2, 0>() == 0.0));
    A_matrix_3x1_fundamentals_type m_A_z_cross_x = cross(m_A_z_dir, m_A_x_dir);
    BOOST_CHECK((m_A_z_cross_x.at<0, 0>() == 0.0));
    BOOST_CHECK((m_A_z_cross_x.at<1, 0>() == 1.0));
    BOOST_CHECK((m_A_z_cross_x.at<2, 0>() == 0.0));
    m_A_z_cross_x = m_A_z_dir ^ m_A_x_dir;
    BOOST_CHECK((m_A_z_cross_x.at<0, 0>() == 0.0));
    BOOST_CHECK((m_A_z_cross_x.at<1, 0>() == 1.0));
    BOOST_CHECK((m_A_z_cross_x.at<2, 0>() == 0.0));
    A_matrix_3x1_fundamentals_type m_A_y_cross_z = cross(m_A_y_dir, m_A_z_dir);
    BOOST_CHECK((m_A_y_cross_z.at<0, 0>() == 1.0));
    BOOST_CHECK((m_A_y_cross_z.at<1, 0>() == 0.0));
    BOOST_CHECK((m_A_y_cross_z.at<2, 0>() == 0.0));
    m_A_y_cross_z = m_A_y_dir ^ m_A_z_dir;
    BOOST_CHECK((m_A_y_cross_z.at<0, 0>() == 1.0));
    BOOST_CHECK((m_A_y_cross_z.at<1, 0>() == 0.0));
    BOOST_CHECK((m_A_y_cross_z.at<2, 0>() == 0.0));
    A_matrix_3x1_fundamentals_type m_A_z_cross_y = cross(m_A_z_dir, m_A_y_dir);
    BOOST_CHECK((m_A_z_cross_y.at<0, 0>() == -1.0));
    BOOST_CHECK((m_A_z_cross_y.at<1, 0>() == 0.0));
    BOOST_CHECK((m_A_z_cross_y.at<2, 0>() == 0.0));
    m_A_z_cross_y = m_A_z_dir ^ m_A_y_dir;
    BOOST_CHECK((m_A_z_cross_y.at<0, 0>() == -1.0));
    BOOST_CHECK((m_A_z_cross_y.at<1, 0>() == 0.0));
    BOOST_CHECK((m_A_z_cross_y.at<2, 0>() == 0.0));

    B_matrix_1x3_fundamentals_type m_B_x_dir;
    m_B_x_dir.at<0, 0>() = 1.0;
    m_B_x_dir.at<0, 1>() = 0.0;
    m_B_x_dir.at<0, 2>() = 0.0;
    B_matrix_1x3_fundamentals_type m_B_y_dir;
    m_B_y_dir.at<0, 0>() = 0.0;
    m_B_y_dir.at<0, 1>() = 1.0;
    m_B_y_dir.at<0, 2>() = 0.0;
    B_matrix_1x3_fundamentals_type m_B_z_dir;
    m_B_z_dir.at<0, 0>() = 0.0;
    m_B_z_dir.at<0, 1>() = 0.0;
    m_B_z_dir.at<0, 2>() = 1.0;

    B_matrix_1x3_fundamentals_type m_B_x_cross_y = cross(m_B_x_dir, m_B_y_dir);
    BOOST_CHECK((m_B_x_cross_y.at<0, 0>() == 0.0));
    BOOST_CHECK((m_B_x_cross_y.at<0, 1>() == 0.0));
    BOOST_CHECK((m_B_x_cross_y.at<0, 2>() == 1.0));
    m_B_x_cross_y = m_B_x_dir ^ m_B_y_dir;
    BOOST_CHECK((m_B_x_cross_y.at<0, 0>() == 0.0));
    BOOST_CHECK((m_B_x_cross_y.at<0, 1>() == 0.0));
    BOOST_CHECK((m_B_x_cross_y.at<0, 2>() == 1.0));
    B_matrix_1x3_fundamentals_type m_B_y_cross_x = cross(m_B_y_dir, m_B_x_dir);
    BOOST_CHECK((m_B_y_cross_x.at<0, 0>() == 0.0));
    BOOST_CHECK((m_B_y_cross_x.at<0, 1>() == 0.0));
    BOOST_CHECK((m_B_y_cross_x.at<0, 2>() == -1.0));
    m_B_y_cross_x = m_B_y_dir ^ m_B_x_dir;
    BOOST_CHECK((m_B_y_cross_x.at<0, 0>() == 0.0));
    BOOST_CHECK((m_B_y_cross_x.at<0, 1>() == 0.0));
    BOOST_CHECK((m_B_y_cross_x.at<0, 2>() == -1.0));
    B_matrix_1x3_fundamentals_type m_B_x_cross_z = cross(m_B_x_dir, m_B_z_dir);
    BOOST_CHECK((m_B_x_cross_z.at<0, 0>() == 0.0));
    BOOST_CHECK((m_B_x_cross_z.at<0, 1>() == -1.0));
    BOOST_CHECK((m_B_x_cross_z.at<0, 2>() == 0.0));
    m_B_x_cross_z = m_B_x_dir ^ m_B_z_dir;
    BOOST_CHECK((m_B_x_cross_z.at<0, 0>() == 0.0));
    BOOST_CHECK((m_B_x_cross_z.at<0, 1>() == -1.0));
    BOOST_CHECK((m_B_x_cross_z.at<0, 2>() == 0.0));
    B_matrix_1x3_fundamentals_type m_B_z_cross_x = cross(m_B_z_dir, m_B_x_dir);
    BOOST_CHECK((m_B_z_cross_x.at<0, 0>() == 0.0));
    BOOST_CHECK((m_B_z_cross_x.at<0, 1>() == 1.0));
    BOOST_CHECK((m_B_z_cross_x.at<0, 2>() == 0.0));
    m_B_z_cross_x = m_B_z_dir ^ m_B_x_dir;
    BOOST_CHECK((m_B_z_cross_x.at<0, 0>() == 0.0));
    BOOST_CHECK((m_B_z_cross_x.at<0, 1>() == 1.0));
    BOOST_CHECK((m_B_z_cross_x.at<0, 2>() == 0.0));
    B_matrix_1x3_fundamentals_type m_B_y_cross_z = cross(m_B_y_dir, m_B_z_dir);
    BOOST_CHECK((m_B_y_cross_z.at<0, 0>() == 1.0));
    BOOST_CHECK((m_B_y_cross_z.at<0, 1>() == 0.0));
    BOOST_CHECK((m_B_y_cross_z.at<0, 2>() == 0.0));
    m_B_y_cross_z = m_B_y_dir ^ m_B_z_dir;
    BOOST_CHECK((m_B_y_cross_z.at<0, 0>() == 1.0));
    BOOST_CHECK((m_B_y_cross_z.at<0, 1>() == 0.0));
    BOOST_CHECK((m_B_y_cross_z.at<0, 2>() == 0.0));
    B_matrix_1x3_fundamentals_type m_B_z_cross_y = cross(m_B_z_dir, m_B_y_dir);
    BOOST_CHECK((m_B_z_cross_y.at<0, 0>() == -1.0));
    BOOST_CHECK((m_B_z_cross_y.at<0, 1>() == 0.0));
    BOOST_CHECK((m_B_z_cross_y.at<0, 2>() == 0.0));
    m_B_z_cross_y = m_B_z_dir ^ m_B_y_dir;
    BOOST_CHECK((m_B_z_cross_y.at<0, 0>() == -1.0));
    BOOST_CHECK((m_B_z_cross_y.at<0, 1>() == 0.0));
    BOOST_CHECK((m_B_z_cross_y.at<0, 2>() == 0.0));

    A_matrix_3x1_fundamentals_type m_A3;
    m_A3.at<0, 0>() = 0.0;
    m_A3.at<1, 0>() = -1.0;
    m_A3.at<2, 0>() = -1.0;
    BOOST_CHECK((sum(m_A3) == -2.0));
    BOOST_CHECK((norm_1(m_A3) == 2.0));
    BOOST_CHECK((norm_2(m_A3) == std::sqrt(2.0)));
    BOOST_CHECK((norm_inf(m_A3) == 1.0));
    BOOST_CHECK((norm_inf_index(m_A3) == 1));

    B_matrix_1x3_fundamentals_type m_B3;
    m_B3.at<0, 0>() = 0.0;
    m_B3.at<0, 1>() = -1.0;
    m_B3.at<0, 2>() = -1.0;
    BOOST_CHECK((sum(m_B3) == -2.0));
    BOOST_CHECK((norm_1(m_B3) == 2.0));
    BOOST_CHECK((norm_2(m_B3) == std::sqrt(2.0)));
    BOOST_CHECK((norm_inf(m_B3) == 1.0));
    BOOST_CHECK((norm_inf_index(m_B3) == 1));

    E_matrix_1x1_fundamentals_type m_E;
    m_E.at<0, 0>() = -2.0;
    BOOST_CHECK((sum(m_E) == -2.0));
    BOOST_CHECK((norm_1(m_E) == 2.0));
    BOOST_CHECK((norm_2(m_E) == std::sqrt(4.0)));
    BOOST_CHECK((norm_inf(m_E) == 2.0));
    BOOST_CHECK((norm_inf_index(m_E) == 0));


    // unit types

    A_matrix_3x1_units_type m_A1_u;
    m_A1_u.at<0, 0>() = time_::from_value(1.0);
    m_A1_u.at<1, 0>() = length::from_value(2.0);
    m_A1_u.at<2, 0>() = dimensionless::from_value(3.0);
    A_matrix_3x1_units_type_2 m_A1_u_2;
    m_A1_u_2.at<0, 0>() = length::from_value(1.0);
    m_A1_u_2.at<1, 0>() = time_::from_value(2.0);
    m_A1_u_2.at<2, 0>() = time_length::from_value(3.0);
    BOOST_CHECK((dot(m_A1_u, m_A1_u_2) == time_length::from_value(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0)));
    BOOST_CHECK((m_A1_u * m_A1_u_2 == time_length::from_value(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0)));

    B_matrix_1x3_units_type m_B1_u;
    m_B1_u.at<0, 0>() = time_::from_value(1.0);
    m_B1_u.at<0, 1>() = time_::from_value(2.0);
    m_B1_u.at<0, 2>() = time_::from_value(3.0);
    B_matrix_1x3_units_type m_B2_u = m_B1_u;
    BOOST_CHECK((dot(m_B1_u, m_B2_u) == time_sq::from_value(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0)));
    BOOST_CHECK((m_B1_u * m_B2_u == time_sq::from_value(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0)));

    A_matrix_3x1_units_type m_A_x_dir_u;
    m_A_x_dir_u.at<0, 0>() = time_::from_value(1.0);
    m_A_x_dir_u.at<1, 0>() = length::from_value(0.0);
    m_A_x_dir_u.at<2, 0>() = dimensionless::from_value(0.0);
    A_matrix_3x1_units_type m_A_y_dir_u;
    m_A_y_dir_u.at<0, 0>() = time_::from_value(0.0);
    m_A_y_dir_u.at<1, 0>() = length::from_value(1.0);
    m_A_y_dir_u.at<2, 0>() = dimensionless::from_value(0.0);
    A_matrix_3x1_units_type m_A_z_dir_u;
    m_A_z_dir_u.at<0, 0>() = time_::from_value(0.0);
    m_A_z_dir_u.at<1, 0>() = length::from_value(0.0);
    m_A_z_dir_u.at<2, 0>() = dimensionless::from_value(1.0);

    typedef bub::result_of::cross_product<
        A_matrix_3x1_units_type,
        A_matrix_3x1_units_type
    >::type A_units_cross_type;
    A_units_cross_type m_A_x_cross_y_u = cross(m_A_x_dir_u, m_A_y_dir_u);
    BOOST_CHECK((m_A_x_cross_y_u.at<0, 0>() == length::from_value(0.0)));
    BOOST_CHECK((m_A_x_cross_y_u.at<1, 0>() == time_::from_value(0.0)));
    BOOST_CHECK((m_A_x_cross_y_u.at<2, 0>() == time_length::from_value(1.0)));
    m_A_x_cross_y_u = m_A_x_dir_u ^ m_A_y_dir_u;
    BOOST_CHECK((m_A_x_cross_y_u.at<0, 0>() == length::from_value(0.0)));
    BOOST_CHECK((m_A_x_cross_y_u.at<1, 0>() == time_::from_value(0.0)));
    BOOST_CHECK((m_A_x_cross_y_u.at<2, 0>() == time_length::from_value(1.0)));
    A_units_cross_type m_A_y_cross_x_u = cross(m_A_y_dir_u, m_A_x_dir_u);
    BOOST_CHECK((m_A_y_cross_x_u.at<0, 0>() == length::from_value(0.0)));
    BOOST_CHECK((m_A_y_cross_x_u.at<1, 0>() == time_::from_value(0.0)));
    BOOST_CHECK((m_A_y_cross_x_u.at<2, 0>() == time_length::from_value(-1.0)));
    m_A_y_cross_x_u = m_A_y_dir_u ^ m_A_x_dir_u;
    BOOST_CHECK((m_A_y_cross_x_u.at<0, 0>() == length::from_value(0.0)));
    BOOST_CHECK((m_A_y_cross_x_u.at<1, 0>() == time_::from_value(0.0)));
    BOOST_CHECK((m_A_y_cross_x_u.at<2, 0>() == time_length::from_value(-1.0)));
    A_units_cross_type m_A_x_cross_z_u = cross(m_A_x_dir_u, m_A_z_dir_u);
    BOOST_CHECK((m_A_x_cross_z_u.at<0, 0>() == length::from_value(0.0)));
    BOOST_CHECK((m_A_x_cross_z_u.at<1, 0>() == time_::from_value(-1.0)));
    BOOST_CHECK((m_A_x_cross_z_u.at<2, 0>() == time_length::from_value(0.0)));
    m_A_x_cross_z_u = m_A_x_dir_u ^ m_A_z_dir_u;
    BOOST_CHECK((m_A_x_cross_z_u.at<0, 0>() == length::from_value(0.0)));
    BOOST_CHECK((m_A_x_cross_z_u.at<1, 0>() == time_::from_value(-1.0)));
    BOOST_CHECK((m_A_x_cross_z_u.at<2, 0>() == time_length::from_value(0.0)));
    A_units_cross_type m_A_z_cross_x_u = cross(m_A_z_dir_u, m_A_x_dir_u);
    BOOST_CHECK((m_A_z_cross_x_u.at<0, 0>() == length::from_value(0.0)));
    BOOST_CHECK((m_A_z_cross_x_u.at<1, 0>() == time_::from_value(1.0)));
    BOOST_CHECK((m_A_z_cross_x_u.at<2, 0>() == time_length::from_value(0.0)));
    m_A_z_cross_x_u = m_A_z_dir_u ^ m_A_x_dir_u;
    BOOST_CHECK((m_A_z_cross_x_u.at<0, 0>() == length::from_value(0.0)));
    BOOST_CHECK((m_A_z_cross_x_u.at<1, 0>() == time_::from_value(1.0)));
    BOOST_CHECK((m_A_z_cross_x_u.at<2, 0>() == time_length::from_value(0.0)));
    A_units_cross_type m_A_y_cross_z_u = cross(m_A_y_dir_u, m_A_z_dir_u);
    BOOST_CHECK((m_A_y_cross_z_u.at<0, 0>() == length::from_value(1.0)));
    BOOST_CHECK((m_A_y_cross_z_u.at<1, 0>() == time_::from_value(0.0)));
    BOOST_CHECK((m_A_y_cross_z_u.at<2, 0>() == time_length::from_value(0.0)));
    m_A_y_cross_z_u = m_A_y_dir_u ^ m_A_z_dir_u;
    BOOST_CHECK((m_A_y_cross_z_u.at<0, 0>() == length::from_value(1.0)));
    BOOST_CHECK((m_A_y_cross_z_u.at<1, 0>() == time_::from_value(0.0)));
    BOOST_CHECK((m_A_y_cross_z_u.at<2, 0>() == time_length::from_value(0.0)));
    A_units_cross_type m_A_z_cross_y_u = cross(m_A_z_dir_u, m_A_y_dir_u);
    BOOST_CHECK((m_A_z_cross_y_u.at<0, 0>() == length::from_value(-1.0)));
    BOOST_CHECK((m_A_z_cross_y_u.at<1, 0>() == time_::from_value(0.0)));
    BOOST_CHECK((m_A_z_cross_y_u.at<2, 0>() == time_length::from_value(0.0)));
    m_A_z_cross_y_u = m_A_z_dir_u ^ m_A_y_dir_u;
    BOOST_CHECK((m_A_z_cross_y_u.at<0, 0>() == length::from_value(-1.0)));
    BOOST_CHECK((m_A_z_cross_y_u.at<1, 0>() == time_::from_value(0.0)));
    BOOST_CHECK((m_A_z_cross_y_u.at<2, 0>() == time_length::from_value(0.0)));

    B_matrix_1x3_units_type m_B_x_dir_u;
    m_B_x_dir_u.at<0, 0>() = time_::from_value(1.0);
    m_B_x_dir_u.at<0, 1>() = time_::from_value(0.0);
    m_B_x_dir_u.at<0, 2>() = time_::from_value(0.0);
    B_matrix_1x3_units_type m_B_y_dir_u;
    m_B_y_dir_u.at<0, 0>() = time_::from_value(0.0);
    m_B_y_dir_u.at<0, 1>() = time_::from_value(1.0);
    m_B_y_dir_u.at<0, 2>() = time_::from_value(0.0);
    B_matrix_1x3_units_type m_B_z_dir_u;
    m_B_z_dir_u.at<0, 0>() = time_::from_value(0.0);
    m_B_z_dir_u.at<0, 1>() = time_::from_value(0.0);
    m_B_z_dir_u.at<0, 2>() = time_::from_value(1.0);

    typedef bub::result_of::cross_product<
        B_matrix_1x3_units_type,
        B_matrix_1x3_units_type
    >::type B_units_cross_type;
    B_units_cross_type m_B_x_cross_y_u = cross(m_B_x_dir_u, m_B_y_dir_u);
    BOOST_CHECK((m_B_x_cross_y_u.at<0, 0>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_x_cross_y_u.at<0, 1>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_x_cross_y_u.at<0, 2>() == time_sq::from_value(1.0)));
    m_B_x_cross_y_u = m_B_x_dir_u ^ m_B_y_dir_u;
    BOOST_CHECK((m_B_x_cross_y_u.at<0, 0>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_x_cross_y_u.at<0, 1>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_x_cross_y_u.at<0, 2>() == time_sq::from_value(1.0)));
    B_units_cross_type m_B_y_cross_x_u = cross(m_B_y_dir_u, m_B_x_dir_u);
    BOOST_CHECK((m_B_y_cross_x_u.at<0, 0>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_y_cross_x_u.at<0, 1>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_y_cross_x_u.at<0, 2>() == time_sq::from_value(-1.0)));
    m_B_y_cross_x_u = m_B_y_dir_u ^ m_B_x_dir_u;
    BOOST_CHECK((m_B_y_cross_x_u.at<0, 0>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_y_cross_x_u.at<0, 1>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_y_cross_x_u.at<0, 2>() == time_sq::from_value(-1.0)));
    B_units_cross_type m_B_x_cross_z_u = cross(m_B_x_dir_u, m_B_z_dir_u);
    BOOST_CHECK((m_B_x_cross_z_u.at<0, 0>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_x_cross_z_u.at<0, 1>() == time_sq::from_value(-1.0)));
    BOOST_CHECK((m_B_x_cross_z_u.at<0, 2>() == time_sq::from_value(0.0)));
    m_B_x_cross_z_u = m_B_x_dir_u ^ m_B_z_dir_u;
    BOOST_CHECK((m_B_x_cross_z_u.at<0, 0>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_x_cross_z_u.at<0, 1>() == time_sq::from_value(-1.0)));
    BOOST_CHECK((m_B_x_cross_z_u.at<0, 2>() == time_sq::from_value(0.0)));
    B_units_cross_type m_B_z_cross_x_u = cross(m_B_z_dir_u, m_B_x_dir_u);
    BOOST_CHECK((m_B_z_cross_x_u.at<0, 0>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_z_cross_x_u.at<0, 1>() == time_sq::from_value(1.0)));
    BOOST_CHECK((m_B_z_cross_x_u.at<0, 2>() == time_sq::from_value(0.0)));
    m_B_z_cross_x_u = m_B_z_dir_u ^ m_B_x_dir_u;
    BOOST_CHECK((m_B_z_cross_x_u.at<0, 0>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_z_cross_x_u.at<0, 1>() == time_sq::from_value(1.0)));
    BOOST_CHECK((m_B_z_cross_x_u.at<0, 2>() == time_sq::from_value(0.0)));
    B_units_cross_type m_B_y_cross_z_u = cross(m_B_y_dir_u, m_B_z_dir_u);
    BOOST_CHECK((m_B_y_cross_z_u.at<0, 0>() == time_sq::from_value(1.0)));
    BOOST_CHECK((m_B_y_cross_z_u.at<0, 1>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_y_cross_z_u.at<0, 2>() == time_sq::from_value(0.0)));
    m_B_y_cross_z_u = m_B_y_dir_u ^ m_B_z_dir_u;
    BOOST_CHECK((m_B_y_cross_z_u.at<0, 0>() == time_sq::from_value(1.0)));
    BOOST_CHECK((m_B_y_cross_z_u.at<0, 1>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_y_cross_z_u.at<0, 2>() == time_sq::from_value(0.0)));
    B_units_cross_type m_B_z_cross_y_u = cross(m_B_z_dir_u, m_B_y_dir_u);
    BOOST_CHECK((m_B_z_cross_y_u.at<0, 0>() == time_sq::from_value(-1.0)));
    BOOST_CHECK((m_B_z_cross_y_u.at<0, 1>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_z_cross_y_u.at<0, 2>() == time_sq::from_value(0.0)));
    m_B_z_cross_y_u = m_B_z_dir_u ^ m_B_y_dir_u;
    BOOST_CHECK((m_B_z_cross_y_u.at<0, 0>() == time_sq::from_value(-1.0)));
    BOOST_CHECK((m_B_z_cross_y_u.at<0, 1>() == time_sq::from_value(0.0)));
    BOOST_CHECK((m_B_z_cross_y_u.at<0, 2>() == time_sq::from_value(0.0)));

    C_matrix_3x1_units_type m_C1_u;
    m_C1_u.at<0, 0>() = dimensionless::from_value(0.0);
    m_C1_u.at<1, 0>() = dimensionless::from_value(-1.0);
    m_C1_u.at<2, 0>() = dimensionless::from_value(-1.0);
    BOOST_CHECK((sum(m_C1_u) == dimensionless::from_value(-2.0)));
    BOOST_CHECK((norm_1(m_C1_u) == dimensionless::from_value(2.0)));
    BOOST_CHECK((norm_2(m_C1_u) == dimensionless::from_value(std::sqrt(2.0))));
    BOOST_CHECK((norm_inf(m_C1_u) == dimensionless::from_value(1.0)));
    BOOST_CHECK((norm_inf_index(m_C1_u) == dimensionless::from_value(1)));

    B_matrix_1x3_units_type m_B3_u;
    m_B3_u.at<0, 0>() = time_::from_value(0.0);
    m_B3_u.at<0, 1>() = time_::from_value(-1.0);
    m_B3_u.at<0, 2>() = time_::from_value(-1.0);
    BOOST_CHECK((sum(m_B3_u) == time_::from_value(-2.0)));
    BOOST_CHECK((norm_1(m_B3_u) == time_::from_value(2.0)));
    BOOST_CHECK((norm_2(m_B3_u) == time_::from_value(std::sqrt(2.0))));
    BOOST_CHECK((norm_inf(m_B3_u) == time_::from_value(1.0)));
    BOOST_CHECK((norm_inf_index(m_B3_u) == 1));

    E_matrix_1x1_units_type m_E_u;
    m_E_u.at<0, 0>() = length::from_value(-2.0);
    BOOST_CHECK((sum(m_E_u) == length::from_value(-2.0)));
    BOOST_CHECK((norm_1(m_E_u) == length::from_value(2.0)));
    BOOST_CHECK((norm_2(m_E_u) == length::from_value(std::sqrt(4.0))));
    BOOST_CHECK((norm_inf(m_E_u) == length::from_value(2.0)));
    BOOST_CHECK((norm_inf_index(m_E_u) == 0));

    // unit types and don't-cares

    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<time_>,
            boost::fusion::vector<bub::_>,
            boost::fusion::vector<dimensionless>
        >
    >::type matrix_3x1_with_dont_care_y;
    matrix_3x1_with_dont_care_y m_A1_u_dc_y;
    m_A1_u_dc_y.at<0, 0>() = time_::from_value(1.0);
    m_A1_u_dc_y.at<2, 0>() = dimensionless::from_value(3.0);

    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<length>,
            boost::fusion::vector<bub::_>,
            boost::fusion::vector<time_length>
        >
    >::type matrix_3x1_with_dont_care_y_2;
    matrix_3x1_with_dont_care_y_2 m_A1_u_dc_y_2;
    m_A1_u_dc_y_2.at<0, 0>() = length::from_value(1.0);
    m_A1_u_dc_y_2.at<2, 0>() = time_length::from_value(3.0);
    BOOST_CHECK((dot(m_A1_u_dc_y, m_A1_u_dc_y_2) == time_length::from_value(1.0 * 1.0 + 3.0 * 3.0)));
    BOOST_CHECK((m_A1_u_dc_y * m_A1_u_dc_y_2 == time_length::from_value(1.0 * 1.0 + 3.0 * 3.0)));

    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<bub::_, time_, time_>
        >
    >::type matrix_1x3_units_type_with_dont_care_x;
    matrix_1x3_units_type_with_dont_care_x m_B_y_dir_u_dc_x;
    m_B_y_dir_u_dc_x.at<0, 1>() = time_::from_value(1.0);
    m_B_y_dir_u_dc_x.at<0, 2>() = time_::from_value(0.0);
    matrix_1x3_units_type_with_dont_care_x m_B_z_dir_u_dc_x;
    m_B_z_dir_u_dc_x.at<0, 1>() = time_::from_value(0.0);
    m_B_z_dir_u_dc_x.at<0, 2>() = time_::from_value(1.0);

    // this cross product will be all don't-cares
    typedef bub::result_of::cross_product<
        matrix_1x3_units_type_with_dont_care_x,
        matrix_1x3_units_type_with_dont_care_x
    >::type matrix_1x3_with_dont_care_x_cross_type;
    matrix_1x3_with_dont_care_x_cross_type m_B_y_cross_z_u_dc_x = cross(m_B_y_dir_u_dc_x, m_B_z_dir_u_dc_x);
    m_B_y_cross_z_u_dc_x = m_B_y_dir_u_dc_x ^ m_B_z_dir_u_dc_x;
    matrix_1x3_with_dont_care_x_cross_type m_B_z_cross_y_u_dc_x = cross(m_B_z_dir_u_dc_x, m_B_y_dir_u_dc_x);
    m_B_z_cross_y_u_dc_x = m_B_z_dir_u_dc_x ^ m_B_y_dir_u_dc_x;

    typedef bub::make_matrix<
        boost::fusion::vector<
            boost::fusion::vector<dimensionless>,
            boost::fusion::vector<dimensionless>,
            boost::fusion::vector<bub::_>
        >
    >::type matrix_3x1_with_dont_care_z;
    matrix_3x1_with_dont_care_z m_C1_u_dc_z;
    m_C1_u_dc_z.at<0, 0>() = dimensionless::from_value(-1.0);
    m_C1_u_dc_z.at<1, 0>() = dimensionless::from_value(-1.0);
    BOOST_CHECK((sum(m_C1_u) == dimensionless::from_value(-2.0)));
    BOOST_CHECK((norm_1(m_C1_u) == dimensionless::from_value(2.0)));
    BOOST_CHECK((norm_2(m_C1_u) == dimensionless::from_value(std::sqrt(2.0))));
    BOOST_CHECK((norm_inf(m_C1_u) == dimensionless::from_value(1.0)));
    BOOST_CHECK((norm_inf_index(m_C1_u) == dimensionless::from_value(1)));

    return 0;
}
