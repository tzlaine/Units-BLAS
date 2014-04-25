// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "operations_tests.hpp"

#include <boost/test/test_tools.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>


typedef bub::matrix<
    std::tuple<double>,
    std::tuple<double>,
    std::tuple<double>,
    std::tuple<double>
> x_vector_type;

struct derived_from_x_vector_type :
    x_vector_type
{};

struct derived_from_D_matrix_4x4_fundamentals_type :
    D_matrix_4x4_fundamentals_type
{};

using b_fundamentals_vector_type = decltype(
    bub::prod(std::declval<D_matrix_4x4_fundamentals_type>(),
              std::declval<x_vector_type>())
);

struct derived_from_b_fundamentals_vector_type :
    b_fundamentals_vector_type
{};

struct derived_from_D_matrix_4x4_units_type :
    D_matrix_4x4_units_type
{};

using b_units_vector_type = decltype(
    bub::prod(std::declval<D_matrix_4x4_units_type>(),
              std::declval<x_vector_type>())
);

struct derived_from_b_units_vector_type :
    b_units_vector_type
{};

int test_main (int, char *[])
{
    double const epsilon = 1.0e-12;

    D_matrix_4x4_fundamentals_type m_D;
    m_D.at<0, 0>() = -1.0;
    m_D.at<0, 1>() = 2.0;
    m_D.at<0, 2>() = 3.0;
    m_D.at<0, 3>() = 4.0;
    m_D.at<1, 0>() = 5.0;
    m_D.at<1, 1>() = 6.0;
    m_D.at<1, 2>() = 7.0;
    m_D.at<1, 3>() = 8.0;
    m_D.at<2, 0>() = 9.0;
    m_D.at<2, 1>() = 10.0;
    m_D.at<2, 2>() = 10.0;
    m_D.at<2, 3>() = 12.0;
    m_D.at<3, 0>() = 13.0;
    m_D.at<3, 1>() = 14.0;
    m_D.at<3, 2>() = 15.0;
    m_D.at<3, 3>() = 16.0;

    derived_from_D_matrix_4x4_fundamentals_type m_D_d;
    m_D_d.at<0, 0>() = -1.0;
    m_D_d.at<0, 1>() = 2.0;
    m_D_d.at<0, 2>() = 3.0;
    m_D_d.at<0, 3>() = 4.0;
    m_D_d.at<1, 0>() = 5.0;
    m_D_d.at<1, 1>() = 6.0;
    m_D_d.at<1, 2>() = 7.0;
    m_D_d.at<1, 3>() = 8.0;
    m_D_d.at<2, 0>() = 9.0;
    m_D_d.at<2, 1>() = 10.0;
    m_D_d.at<2, 2>() = 10.0;
    m_D_d.at<2, 3>() = 12.0;
    m_D_d.at<3, 0>() = 13.0;
    m_D_d.at<3, 1>() = 14.0;
    m_D_d.at<3, 2>() = 15.0;
    m_D_d.at<3, 3>() = 16.0;

    BOOST_CHECK_CLOSE(determinant(m_D), -32.0, epsilon);

    BOOST_CHECK_CLOSE(determinant(m_D_d), -32.0, epsilon);

#if 0
    typedef bub::result_of::inverse<
        D_matrix_4x4_fundamentals_type
    >::type D_matrix_4x4_fundamentals_type_inverse;
    D_matrix_4x4_fundamentals_type_inverse m_D_inv;
    m_D_inv = inverse(m_D);
    BOOST_CHECK_CLOSE((m_D_inv.at<0, 0>()), -0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<0, 1>()), 0.75, epsilon);
    BOOST_CHECK_SMALL((m_D_inv.at<0, 2>()), epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<0, 3>()), -0.25, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<1, 0>()), 0.75, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<1, 1>()), -2.375, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<1, 2>()), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<1, 3>()), 0.625, epsilon);
    BOOST_CHECK_SMALL((m_D_inv.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<2, 1>()), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<2, 2>()), -1.0, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<2, 3>()), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<3, 0>()), -0.25, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<3, 1>()), 1.0, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<3, 2>()), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<3, 3>()), -0.75, epsilon);

    m_D_inv = inverse(m_D_d);
    BOOST_CHECK_CLOSE((m_D_inv.at<0, 0>()), -0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<0, 1>()), 0.75, epsilon);
    BOOST_CHECK_SMALL((m_D_inv.at<0, 2>()), epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<0, 3>()), -0.25, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<1, 0>()), 0.75, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<1, 1>()), -2.375, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<1, 2>()), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<1, 3>()), 0.625, epsilon);
    BOOST_CHECK_SMALL((m_D_inv.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<2, 1>()), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<2, 2>()), -1.0, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<2, 3>()), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<3, 0>()), -0.25, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<3, 1>()), 1.0, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<3, 2>()), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_inv.at<3, 3>()), -0.75, epsilon);

    b_fundamentals_vector_type m_b;
    m_b.at<0, 0>() = 1.0;
    m_b.at<1, 0>() = 1.0;
    m_b.at<2, 0>() = 1.0;
    m_b.at<3, 0>() = 1.0;

    derived_from_b_fundamentals_vector_type m_b_d;
    m_b_d.at<0, 0>() = 1.0;
    m_b_d.at<1, 0>() = 1.0;
    m_b_d.at<2, 0>() = 1.0;
    m_b_d.at<3, 0>() = 1.0;

    x_vector_type m_x;
    solve(m_D, m_b, m_x);
    BOOST_CHECK_SMALL((m_x.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x.at<3, 0>()), 0.5, epsilon);

    solve(m_D_d, m_b, m_x);
    BOOST_CHECK_SMALL((m_x.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x.at<3, 0>()), 0.5, epsilon);

    solve(m_D, m_b_d, m_x);
    BOOST_CHECK_SMALL((m_x.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x.at<3, 0>()), 0.5, epsilon);

    solve(m_D_d, m_b_d, m_x);
    BOOST_CHECK_SMALL((m_x.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x.at<3, 0>()), 0.5, epsilon);

    derived_from_x_vector_type m_x_d;
    solve(m_D, m_b, m_x_d);
    BOOST_CHECK_SMALL((m_x_d.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_d.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_d.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_d.at<3, 0>()), 0.5, epsilon);

    solve(m_D_d, m_b, m_x_d);
    BOOST_CHECK_SMALL((m_x_d.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_d.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_d.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_d.at<3, 0>()), 0.5, epsilon);

    solve(m_D, m_b_d, m_x_d);
    BOOST_CHECK_SMALL((m_x_d.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_d.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_d.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_d.at<3, 0>()), 0.5, epsilon);

    solve(m_D_d, m_b_d, m_x_d);
    BOOST_CHECK_SMALL((m_x_d.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_d.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_d.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_d.at<3, 0>()), 0.5, epsilon);

    D_matrix_4x4_fundamentals_type m_D_zero;
    BOOST_CHECK((determinant(m_D_zero) == 0.0));
    BOOST_CHECK_THROW(inverse(m_D_zero), bub::singular_matrix);
    BOOST_CHECK_THROW(solve(m_D_zero, m_b, m_x), bub::singular_matrix);

    derived_from_D_matrix_4x4_fundamentals_type m_D_zero_d;
    BOOST_CHECK((determinant(m_D_zero_d) == 0.0));
    BOOST_CHECK_THROW(inverse(m_D_zero_d), bub::singular_matrix);
    BOOST_CHECK_THROW(solve(m_D_zero_d, m_b, m_x), bub::singular_matrix);
#endif

    D_matrix_4x4_units_type m_D_u;
    m_D_u.at<0, 0>() = length::from_value(-1.0);
    m_D_u.at<0, 1>() = length::from_value(2.0);
    m_D_u.at<0, 2>() = length::from_value(3.0);
    m_D_u.at<0, 3>() = length::from_value(4.0);
    m_D_u.at<1, 0>() = length::from_value(5.0);
    m_D_u.at<1, 1>() = length::from_value(6.0);
    m_D_u.at<1, 2>() = length::from_value(7.0);
    m_D_u.at<1, 3>() = length::from_value(8.0);
    m_D_u.at<2, 0>() = length::from_value(9.0);
    m_D_u.at<2, 1>() = length::from_value(10.0);
    m_D_u.at<2, 2>() = length::from_value(10.0);
    m_D_u.at<2, 3>() = length::from_value(12.0);
    m_D_u.at<3, 0>() = length::from_value(13.0);
    m_D_u.at<3, 1>() = length::from_value(14.0);
    m_D_u.at<3, 2>() = length::from_value(15.0);
    m_D_u.at<3, 3>() = length::from_value(16.0);

    derived_from_D_matrix_4x4_units_type m_D_u_d;
    m_D_u_d.at<0, 0>() = length::from_value(-1.0);
    m_D_u_d.at<0, 1>() = length::from_value(2.0);
    m_D_u_d.at<0, 2>() = length::from_value(3.0);
    m_D_u_d.at<0, 3>() = length::from_value(4.0);
    m_D_u_d.at<1, 0>() = length::from_value(5.0);
    m_D_u_d.at<1, 1>() = length::from_value(6.0);
    m_D_u_d.at<1, 2>() = length::from_value(7.0);
    m_D_u_d.at<1, 3>() = length::from_value(8.0);
    m_D_u_d.at<2, 0>() = length::from_value(9.0);
    m_D_u_d.at<2, 1>() = length::from_value(10.0);
    m_D_u_d.at<2, 2>() = length::from_value(10.0);
    m_D_u_d.at<2, 3>() = length::from_value(12.0);
    m_D_u_d.at<3, 0>() = length::from_value(13.0);
    m_D_u_d.at<3, 1>() = length::from_value(14.0);
    m_D_u_d.at<3, 2>() = length::from_value(15.0);
    m_D_u_d.at<3, 3>() = length::from_value(16.0);

    BOOST_CHECK_CLOSE(determinant(m_D_u).value(), -32.0, epsilon);

    BOOST_CHECK_CLOSE(determinant(m_D_u_d).value(), -32.0, epsilon);

#if 0
    typedef bub::result_of::inverse<
        D_matrix_4x4_units_type
    >::type D_matrix_4x4_units_type_inverse;
    D_matrix_4x4_units_type_inverse m_D_u_inv;
    m_D_u_inv = inverse(m_D_u);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<0, 0>()).value(), -0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<0, 1>()).value(), 0.75, epsilon);
    BOOST_CHECK_SMALL((m_D_u_inv.at<0, 2>()).value(), epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<0, 3>()).value(), -0.25, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<1, 0>()).value(), 0.75, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<1, 1>()).value(), -2.375, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<1, 2>()).value(), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<1, 3>()).value(), 0.625, epsilon);
    BOOST_CHECK_SMALL((m_D_u_inv.at<2, 0>()).value(), epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<2, 1>()).value(), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<2, 2>()).value(), -1.0, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<2, 3>()).value(), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<3, 0>()).value(), -0.25, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<3, 1>()).value(), 1.0, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<3, 2>()).value(), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<3, 3>()).value(), -0.75, epsilon);

    m_D_u_inv = inverse(m_D_u_d);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<0, 0>()).value(), -0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<0, 1>()).value(), 0.75, epsilon);
    BOOST_CHECK_SMALL((m_D_u_inv.at<0, 2>()).value(), epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<0, 3>()).value(), -0.25, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<1, 0>()).value(), 0.75, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<1, 1>()).value(), -2.375, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<1, 2>()).value(), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<1, 3>()).value(), 0.625, epsilon);
    BOOST_CHECK_SMALL((m_D_u_inv.at<2, 0>()).value(), epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<2, 1>()).value(), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<2, 2>()).value(), -1.0, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<2, 3>()).value(), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<3, 0>()).value(), -0.25, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<3, 1>()).value(), 1.0, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<3, 2>()).value(), 0.5, epsilon);
    BOOST_CHECK_CLOSE((m_D_u_inv.at<3, 3>()).value(), -0.75, epsilon);

    b_units_vector_type m_b_u;
    m_b_u.at<0, 0>() = length::from_value(1.0);
    m_b_u.at<1, 0>() = length::from_value(1.0);
    m_b_u.at<2, 0>() = length::from_value(1.0);
    m_b_u.at<3, 0>() = length::from_value(1.0);

    derived_from_b_units_vector_type m_b_u_d;
    m_b_u_d.at<0, 0>() = length::from_value(1.0);
    m_b_u_d.at<1, 0>() = length::from_value(1.0);
    m_b_u_d.at<2, 0>() = length::from_value(1.0);
    m_b_u_d.at<3, 0>() = length::from_value(1.0);

    x_vector_type m_x_u;
    solve(m_D_u, m_b_u, m_x_u);
    BOOST_CHECK_SMALL((m_x_u.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_u.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u.at<3, 0>()), 0.5, epsilon);

    solve(m_D_u_d, m_b_u, m_x_u);
    BOOST_CHECK_SMALL((m_x_u.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_u.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u.at<3, 0>()), 0.5, epsilon);

    solve(m_D_u, m_b_u_d, m_x_u);
    BOOST_CHECK_SMALL((m_x_u.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_u.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u.at<3, 0>()), 0.5, epsilon);

    solve(m_D_u_d, m_b_u_d, m_x_u);
    BOOST_CHECK_SMALL((m_x_u.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_u.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u.at<3, 0>()), 0.5, epsilon);

    derived_from_x_vector_type m_x_u_d;
    solve(m_D_u, m_b_u, m_x_u_d);
    BOOST_CHECK_SMALL((m_x_u_d.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u_d.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_u_d.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u_d.at<3, 0>()), 0.5, epsilon);

    solve(m_D_u_d, m_b_u, m_x_u_d);
    BOOST_CHECK_SMALL((m_x_u_d.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u_d.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_u_d.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u_d.at<3, 0>()), 0.5, epsilon);

    solve(m_D_u, m_b_u_d, m_x_u_d);
    BOOST_CHECK_SMALL((m_x_u_d.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u_d.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_u_d.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u_d.at<3, 0>()), 0.5, epsilon);

    solve(m_D_u_d, m_b_u_d, m_x_u_d);
    BOOST_CHECK_SMALL((m_x_u_d.at<0, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u_d.at<1, 0>()), -0.5, epsilon);
    BOOST_CHECK_SMALL((m_x_u_d.at<2, 0>()), epsilon);
    BOOST_CHECK_CLOSE((m_x_u_d.at<3, 0>()), 0.5, epsilon);

    D_matrix_4x4_units_type m_D_u_zero;
    BOOST_CHECK((determinant(m_D_u_zero).value() == 0.0));
    BOOST_CHECK_THROW(inverse(m_D_u_zero), bub::singular_matrix);
    BOOST_CHECK_THROW(solve(m_D_u_zero, m_b_u, m_x_u), bub::singular_matrix);

    derived_from_D_matrix_4x4_units_type m_D_u_zero_d;
    BOOST_CHECK((determinant(m_D_u_zero_d).value() == 0.0));
    BOOST_CHECK_THROW(inverse(m_D_u_zero_d), bub::singular_matrix);
    BOOST_CHECK_THROW(solve(m_D_u_zero_d, m_b_u, m_x_u), bub::singular_matrix);
#endif

    return 0;
}
