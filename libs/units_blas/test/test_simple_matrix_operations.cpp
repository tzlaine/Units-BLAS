// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "operations_tests.hpp"

#include <boost/mpl/vector_c.hpp>

#include <boost/test/minimal.hpp>


typedef boost::mpl::vector_c<std::size_t, 0, 1> rows_type_1;
typedef boost::mpl::vector_c<std::size_t, 0, 1> columns_type_1;

// note that rows_type_2 reorders rows 2 and 3
typedef boost::mpl::vector_c<std::size_t, 0, 3, 2> rows_type_2;
typedef boost::mpl::vector_c<std::size_t, 2> columns_type_2;

typedef boost::mpl::vector_c<std::size_t, 0> rows_type_3;
typedef boost::mpl::vector_c<std::size_t, 1, 3> columns_type_3;

struct derived_from_D_matrix_4x4_fundamentals_type :
    D_matrix_4x4_fundamentals_type
{};

struct derived_from_D_matrix_4x4_units_type :
    D_matrix_4x4_units_type
{};

int test_main (int, char *[])
{
    // fundamental types

    A_matrix_3x1_fundamentals_type m_A1;

    bub::at<boost::mpl::size_t<0>, boost::mpl::size_t<0> >(m_A1) = 3.0;
    bub::at<boost::mpl::size_t<1>, boost::mpl::size_t<0> >(m_A1) = 1.0;
    bub::at<boost::mpl::size_t<2>, boost::mpl::size_t<0> >(m_A1) = 2.0;
    BOOST_CHECK((m_A1.at<0, 0>() == 3.0));
    BOOST_CHECK((m_A1.at<1, 0>() == 1.0));
    BOOST_CHECK((m_A1.at<2, 0>() == 2.0));

    bub::at<boost::mpl::size_t<0>, boost::mpl::size_t<0> >(m_A1) = 2.0;
    bub::at<boost::mpl::size_t<1>, boost::mpl::size_t<0> >(m_A1) = 3.0;
    bub::at<boost::mpl::size_t<2>, boost::mpl::size_t<0> >(m_A1) = 1.0;
    BOOST_CHECK((m_A1.at<0, 0>() == 2.0));
    BOOST_CHECK((m_A1.at<1, 0>() == 3.0));
    BOOST_CHECK((m_A1.at<2, 0>() == 1.0));

    bub::at_c<0, 0>(m_A1) = 1.0;
    bub::at_c<1, 0>(m_A1) = 2.0;
    bub::at_c<2, 0>(m_A1) = 3.0;
    BOOST_CHECK((m_A1.at<0, 0>() == 1.0));
    BOOST_CHECK((m_A1.at<1, 0>() == 2.0));
    BOOST_CHECK((m_A1.at<2, 0>() == 3.0));

    A_matrix_3x1_fundamentals_type m_A2 = prod(3.0, m_A1);
    BOOST_CHECK((m_A2.at<0, 0>() == 3.0));
    BOOST_CHECK((m_A2.at<1, 0>() == 6.0));
    BOOST_CHECK((m_A2.at<2, 0>() == 9.0));
    m_A2 = 3.0 * m_A1;
    BOOST_CHECK((m_A2.at<0, 0>() == 3.0));
    BOOST_CHECK((m_A2.at<1, 0>() == 6.0));
    BOOST_CHECK((m_A2.at<2, 0>() == 9.0));

    A_matrix_3x1_fundamentals_type m_A3 = prod(m_A1, 3.0);
    BOOST_CHECK((m_A3.at<0, 0>() == 3.0));
    BOOST_CHECK((m_A3.at<1, 0>() == 6.0));
    BOOST_CHECK((m_A3.at<2, 0>() == 9.0));
    m_A3 = m_A1 * 3.0;
    BOOST_CHECK((m_A3.at<0, 0>() == 3.0));
    BOOST_CHECK((m_A3.at<1, 0>() == 6.0));
    BOOST_CHECK((m_A3.at<2, 0>() == 9.0));

    A_matrix_3x1_fundamentals_type m_A4 = div(m_A1, 3.0);
    BOOST_CHECK((m_A4.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_A4.at<1, 0>() == 2.0 / 3.0));
    BOOST_CHECK((m_A4.at<2, 0>() == 3.0 / 3.0));
    m_A4 = m_A1 / 3.0;
    BOOST_CHECK((m_A4.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_A4.at<1, 0>() == 2.0 / 3.0));
    BOOST_CHECK((m_A4.at<2, 0>() == 3.0 / 3.0));

    A_matrix_3x1_fundamentals_type m_A5 = neg(m_A1);
    BOOST_CHECK((m_A5.at<0, 0>() == -1.0));
    BOOST_CHECK((m_A5.at<1, 0>() == -2.0));
    BOOST_CHECK((m_A5.at<2, 0>() == -3.0));
    m_A5 = -m_A1;
    BOOST_CHECK((m_A5.at<0, 0>() == -1.0));
    BOOST_CHECK((m_A5.at<1, 0>() == -2.0));
    BOOST_CHECK((m_A5.at<2, 0>() == -3.0));

    D_matrix_4x4_fundamentals_type m_D1;
    m_D1.at<0, 0>() = 1.0;
    m_D1.at<0, 1>() = 2.0;
    m_D1.at<0, 2>() = 3.0;
    m_D1.at<0, 3>() = 4.0;
    m_D1.at<1, 0>() = 5.0;
    m_D1.at<1, 1>() = 6.0;
    m_D1.at<1, 2>() = 7.0;
    m_D1.at<1, 3>() = 8.0;
    m_D1.at<2, 0>() = 9.0;
    m_D1.at<2, 1>() = 10.0;
    m_D1.at<2, 2>() = 11.0;
    m_D1.at<2, 3>() = 12.0;
    m_D1.at<3, 0>() = 13.0;
    m_D1.at<3, 1>() = 14.0;
    m_D1.at<3, 2>() = 15.0;
    m_D1.at<3, 3>() = 16.0;

    derived_from_D_matrix_4x4_fundamentals_type m_D1_d;
    m_D1_d.at<0, 0>() = 1.0;
    m_D1_d.at<0, 1>() = 2.0;
    m_D1_d.at<0, 2>() = 3.0;
    m_D1_d.at<0, 3>() = 4.0;
    m_D1_d.at<1, 0>() = 5.0;
    m_D1_d.at<1, 1>() = 6.0;
    m_D1_d.at<1, 2>() = 7.0;
    m_D1_d.at<1, 3>() = 8.0;
    m_D1_d.at<2, 0>() = 9.0;
    m_D1_d.at<2, 1>() = 10.0;
    m_D1_d.at<2, 2>() = 11.0;
    m_D1_d.at<2, 3>() = 12.0;
    m_D1_d.at<3, 0>() = 13.0;
    m_D1_d.at<3, 1>() = 14.0;
    m_D1_d.at<3, 2>() = 15.0;
    m_D1_d.at<3, 3>() = 16.0;

    D_matrix_4x4_fundamentals_type m_D2;
    m_D2 = transpose(m_D1);
    BOOST_CHECK((m_D2.at<0, 0>() == 1.0));
    BOOST_CHECK((m_D2.at<1, 0>() == 2.0));
    BOOST_CHECK((m_D2.at<2, 0>() == 3.0));
    BOOST_CHECK((m_D2.at<3, 0>() == 4.0));
    BOOST_CHECK((m_D2.at<0, 1>() == 5.0));
    BOOST_CHECK((m_D2.at<1, 1>() == 6.0));
    BOOST_CHECK((m_D2.at<2, 1>() == 7.0));
    BOOST_CHECK((m_D2.at<3, 1>() == 8.0));
    BOOST_CHECK((m_D2.at<0, 2>() == 9.0));
    BOOST_CHECK((m_D2.at<1, 2>() == 10.0));
    BOOST_CHECK((m_D2.at<2, 2>() == 11.0));
    BOOST_CHECK((m_D2.at<3, 2>() == 12.0));
    BOOST_CHECK((m_D2.at<0, 3>() == 13.0));
    BOOST_CHECK((m_D2.at<1, 3>() == 14.0));
    BOOST_CHECK((m_D2.at<2, 3>() == 15.0));
    BOOST_CHECK((m_D2.at<3, 3>() == 16.0));

    m_D2 = transpose(m_D1_d);
    BOOST_CHECK((m_D2.at<0, 0>() == 1.0));
    BOOST_CHECK((m_D2.at<1, 0>() == 2.0));
    BOOST_CHECK((m_D2.at<2, 0>() == 3.0));
    BOOST_CHECK((m_D2.at<3, 0>() == 4.0));
    BOOST_CHECK((m_D2.at<0, 1>() == 5.0));
    BOOST_CHECK((m_D2.at<1, 1>() == 6.0));
    BOOST_CHECK((m_D2.at<2, 1>() == 7.0));
    BOOST_CHECK((m_D2.at<3, 1>() == 8.0));
    BOOST_CHECK((m_D2.at<0, 2>() == 9.0));
    BOOST_CHECK((m_D2.at<1, 2>() == 10.0));
    BOOST_CHECK((m_D2.at<2, 2>() == 11.0));
    BOOST_CHECK((m_D2.at<3, 2>() == 12.0));
    BOOST_CHECK((m_D2.at<0, 3>() == 13.0));
    BOOST_CHECK((m_D2.at<1, 3>() == 14.0));
    BOOST_CHECK((m_D2.at<2, 3>() == 15.0));
    BOOST_CHECK((m_D2.at<3, 3>() == 16.0));

    typedef bub::result_of::slice<
        D_matrix_4x4_fundamentals_type,
        rows_type_1,
        columns_type_1
    >::type D_fundamentals_slice_type_1;
    D_fundamentals_slice_type_1 m_D_slice_1;
    m_D_slice_1 = bub::slice<rows_type_1, columns_type_1>(m_D1);
    BOOST_CHECK((m_D_slice_1.at<0, 0>() == 1.0));
    BOOST_CHECK((m_D_slice_1.at<0, 1>() == 2.0));
    BOOST_CHECK((m_D_slice_1.at<1, 0>() == 5.0));
    BOOST_CHECK((m_D_slice_1.at<1, 1>() == 6.0));

    m_D_slice_1 = bub::slice<rows_type_1, columns_type_1>(m_D1_d);
    BOOST_CHECK((m_D_slice_1.at<0, 0>() == 1.0));
    BOOST_CHECK((m_D_slice_1.at<0, 1>() == 2.0));
    BOOST_CHECK((m_D_slice_1.at<1, 0>() == 5.0));
    BOOST_CHECK((m_D_slice_1.at<1, 1>() == 6.0));

    typedef bub::result_of::slice<
        D_matrix_4x4_fundamentals_type,
        rows_type_2,
        columns_type_2
    >::type D_fundamentals_slice_type_2;
    D_fundamentals_slice_type_2 m_D_slice_2;
    m_D_slice_2 = bub::slice<rows_type_2, columns_type_2>(m_D1);
    BOOST_CHECK((m_D_slice_2.at<0, 0>() == 3.0));
    BOOST_CHECK((m_D_slice_2.at<1, 0>() == 15.0));
    BOOST_CHECK((m_D_slice_2.at<2, 0>() == 11.0));

    m_D_slice_2 = bub::slice<rows_type_2, columns_type_2>(m_D1_d);
    BOOST_CHECK((m_D_slice_2.at<0, 0>() == 3.0));
    BOOST_CHECK((m_D_slice_2.at<1, 0>() == 15.0));
    BOOST_CHECK((m_D_slice_2.at<2, 0>() == 11.0));

    typedef bub::result_of::slice<
        D_matrix_4x4_fundamentals_type,
        rows_type_3,
        columns_type_3
    >::type D_fundamentals_slice_type_3;
    D_fundamentals_slice_type_3 m_D_slice_3;
    m_D_slice_3 = bub::slice<rows_type_3, columns_type_3>(m_D1);
    BOOST_CHECK((m_D_slice_3.at<0, 0>() == 2.0));
    BOOST_CHECK((m_D_slice_3.at<0, 1>() == 4.0));

    m_D_slice_3 = bub::slice<rows_type_3, columns_type_3>(m_D1_d);
    BOOST_CHECK((m_D_slice_3.at<0, 0>() == 2.0));
    BOOST_CHECK((m_D_slice_3.at<0, 1>() == 4.0));


    // unit types

    A_matrix_3x1_units_type m_A1_u;

    bub::at<boost::mpl::size_t<0>, boost::mpl::size_t<0> >(m_A1_u) = time_::from_value(3.0);
    bub::at<boost::mpl::size_t<1>, boost::mpl::size_t<0> >(m_A1_u) = length::from_value(1.0);
    bub::at<boost::mpl::size_t<2>, boost::mpl::size_t<0> >(m_A1_u) = dimensionless::from_value(2.0);
    BOOST_CHECK((m_A1_u.at<0, 0>().value() == 3.0));
    BOOST_CHECK((m_A1_u.at<1, 0>().value() == 1.0));
    BOOST_CHECK((m_A1_u.at<2, 0>().value() == 2.0));

    bub::at<boost::mpl::size_t<0>, boost::mpl::size_t<0> >(m_A1_u) = time_::from_value(2.0);
    bub::at<boost::mpl::size_t<1>, boost::mpl::size_t<0> >(m_A1_u) = length::from_value(3.0);
    bub::at<boost::mpl::size_t<2>, boost::mpl::size_t<0> >(m_A1_u) = dimensionless::from_value(1.0);
    BOOST_CHECK((m_A1_u.at<0, 0>().value() == 2.0));
    BOOST_CHECK((m_A1_u.at<1, 0>().value() == 3.0));
    BOOST_CHECK((m_A1_u.at<2, 0>().value() == 1.0));

    bub::at_c<0, 0>(m_A1_u) = time_::from_value(1.0);
    bub::at_c<1, 0>(m_A1_u) = length::from_value(2.0);
    bub::at_c<2, 0>(m_A1_u) = dimensionless::from_value(3.0);
    BOOST_CHECK((m_A1_u.at<0, 0>().value() == 1.0));
    BOOST_CHECK((m_A1_u.at<1, 0>().value() == 2.0));
    BOOST_CHECK((m_A1_u.at<2, 0>().value() == 3.0));

    A_matrix_3x1_units_type m_A2_u = prod(3.0, m_A1_u);
    BOOST_CHECK((m_A2_u.at<0, 0>().value() == 3.0));
    BOOST_CHECK((m_A2_u.at<1, 0>().value() == 6.0));
    BOOST_CHECK((m_A2_u.at<2, 0>().value() == 9.0));
    m_A2_u = 3.0 * m_A1_u;
    BOOST_CHECK((m_A2_u.at<0, 0>().value() == 3.0));
    BOOST_CHECK((m_A2_u.at<1, 0>().value() == 6.0));
    BOOST_CHECK((m_A2_u.at<2, 0>().value() == 9.0));

    A_matrix_3x1_units_type m_A3_u = prod(m_A1_u, 3.0);
    BOOST_CHECK((m_A3_u.at<0, 0>().value() == 3.0));
    BOOST_CHECK((m_A3_u.at<1, 0>().value() == 6.0));
    BOOST_CHECK((m_A3_u.at<2, 0>().value() == 9.0));
    m_A3_u = m_A1_u * 3.0;
    BOOST_CHECK((m_A3_u.at<0, 0>().value() == 3.0));
    BOOST_CHECK((m_A3_u.at<1, 0>().value() == 6.0));
    BOOST_CHECK((m_A3_u.at<2, 0>().value() == 9.0));

    A_matrix_3x1_units_type m_A4_u = div(m_A1_u, 3.0);
    BOOST_CHECK((m_A4_u.at<0, 0>().value() == 1.0 / 3.0));
    BOOST_CHECK((m_A4_u.at<1, 0>().value() == 2.0 / 3.0));
    BOOST_CHECK((m_A4_u.at<2, 0>().value() == 3.0 / 3.0));
    m_A4_u = m_A1_u / 3.0;
    BOOST_CHECK((m_A4_u.at<0, 0>().value() == 1.0 / 3.0));
    BOOST_CHECK((m_A4_u.at<1, 0>().value() == 2.0 / 3.0));
    BOOST_CHECK((m_A4_u.at<2, 0>().value() == 3.0 / 3.0));

    A_matrix_3x1_units_type m_A5_u = neg(m_A1_u);
    BOOST_CHECK((m_A5_u.at<0, 0>().value() == -1.0));
    BOOST_CHECK((m_A5_u.at<1, 0>().value() == -2.0));
    BOOST_CHECK((m_A5_u.at<2, 0>().value() == -3.0));
    m_A5_u = -m_A1_u;
    BOOST_CHECK((m_A5_u.at<0, 0>().value() == -1.0));
    BOOST_CHECK((m_A5_u.at<1, 0>().value() == -2.0));
    BOOST_CHECK((m_A5_u.at<2, 0>().value() == -3.0));

    D_matrix_4x4_units_type m_D1_u;
    m_D1_u.at<0, 0>() = length::from_value(1.0);
    m_D1_u.at<0, 1>() = length::from_value(2.0);
    m_D1_u.at<0, 2>() = length::from_value(3.0);
    m_D1_u.at<0, 3>() = length::from_value(4.0);
    m_D1_u.at<1, 0>() = length::from_value(5.0);
    m_D1_u.at<1, 1>() = length::from_value(6.0);
    m_D1_u.at<1, 2>() = length::from_value(7.0);
    m_D1_u.at<1, 3>() = length::from_value(8.0);
    m_D1_u.at<2, 0>() = length::from_value(9.0);
    m_D1_u.at<2, 1>() = length::from_value(10.0);
    m_D1_u.at<2, 2>() = length::from_value(11.0);
    m_D1_u.at<2, 3>() = length::from_value(12.0);
    m_D1_u.at<3, 0>() = length::from_value(13.0);
    m_D1_u.at<3, 1>() = length::from_value(14.0);
    m_D1_u.at<3, 2>() = length::from_value(15.0);
    m_D1_u.at<3, 3>() = length::from_value(16.0);

    derived_from_D_matrix_4x4_units_type m_D1_u_d;
    m_D1_u_d.at<0, 0>() = length::from_value(1.0);
    m_D1_u_d.at<0, 1>() = length::from_value(2.0);
    m_D1_u_d.at<0, 2>() = length::from_value(3.0);
    m_D1_u_d.at<0, 3>() = length::from_value(4.0);
    m_D1_u_d.at<1, 0>() = length::from_value(5.0);
    m_D1_u_d.at<1, 1>() = length::from_value(6.0);
    m_D1_u_d.at<1, 2>() = length::from_value(7.0);
    m_D1_u_d.at<1, 3>() = length::from_value(8.0);
    m_D1_u_d.at<2, 0>() = length::from_value(9.0);
    m_D1_u_d.at<2, 1>() = length::from_value(10.0);
    m_D1_u_d.at<2, 2>() = length::from_value(11.0);
    m_D1_u_d.at<2, 3>() = length::from_value(12.0);
    m_D1_u_d.at<3, 0>() = length::from_value(13.0);
    m_D1_u_d.at<3, 1>() = length::from_value(14.0);
    m_D1_u_d.at<3, 2>() = length::from_value(15.0);
    m_D1_u_d.at<3, 3>() = length::from_value(16.0);

    D_matrix_4x4_units_type m_D2_u;
    m_D2_u = transpose(m_D1_u);
    BOOST_CHECK((m_D2_u.at<0, 0>().value() == 1.0));
    BOOST_CHECK((m_D2_u.at<1, 0>().value() == 2.0));
    BOOST_CHECK((m_D2_u.at<2, 0>().value() == 3.0));
    BOOST_CHECK((m_D2_u.at<3, 0>().value() == 4.0));
    BOOST_CHECK((m_D2_u.at<0, 1>().value() == 5.0));
    BOOST_CHECK((m_D2_u.at<1, 1>().value() == 6.0));
    BOOST_CHECK((m_D2_u.at<2, 1>().value() == 7.0));
    BOOST_CHECK((m_D2_u.at<3, 1>().value() == 8.0));
    BOOST_CHECK((m_D2_u.at<0, 2>().value() == 9.0));
    BOOST_CHECK((m_D2_u.at<1, 2>().value() == 10.0));
    BOOST_CHECK((m_D2_u.at<2, 2>().value() == 11.0));
    BOOST_CHECK((m_D2_u.at<3, 2>().value() == 12.0));
    BOOST_CHECK((m_D2_u.at<0, 3>().value() == 13.0));
    BOOST_CHECK((m_D2_u.at<1, 3>().value() == 14.0));
    BOOST_CHECK((m_D2_u.at<2, 3>().value() == 15.0));
    BOOST_CHECK((m_D2_u.at<3, 3>().value() == 16.0));

    m_D2_u = transpose(m_D1_u_d);
    BOOST_CHECK((m_D2_u.at<0, 0>().value() == 1.0));
    BOOST_CHECK((m_D2_u.at<1, 0>().value() == 2.0));
    BOOST_CHECK((m_D2_u.at<2, 0>().value() == 3.0));
    BOOST_CHECK((m_D2_u.at<3, 0>().value() == 4.0));
    BOOST_CHECK((m_D2_u.at<0, 1>().value() == 5.0));
    BOOST_CHECK((m_D2_u.at<1, 1>().value() == 6.0));
    BOOST_CHECK((m_D2_u.at<2, 1>().value() == 7.0));
    BOOST_CHECK((m_D2_u.at<3, 1>().value() == 8.0));
    BOOST_CHECK((m_D2_u.at<0, 2>().value() == 9.0));
    BOOST_CHECK((m_D2_u.at<1, 2>().value() == 10.0));
    BOOST_CHECK((m_D2_u.at<2, 2>().value() == 11.0));
    BOOST_CHECK((m_D2_u.at<3, 2>().value() == 12.0));
    BOOST_CHECK((m_D2_u.at<0, 3>().value() == 13.0));
    BOOST_CHECK((m_D2_u.at<1, 3>().value() == 14.0));
    BOOST_CHECK((m_D2_u.at<2, 3>().value() == 15.0));
    BOOST_CHECK((m_D2_u.at<3, 3>().value() == 16.0));

    typedef bub::result_of::slice<
        D_matrix_4x4_units_type,
        rows_type_1,
        columns_type_1
    >::type D_units_slice_type_1;
    D_units_slice_type_1 m_D_slice_1_u;
    m_D_slice_1_u = bub::slice<rows_type_1, columns_type_1>(m_D1_u);
    BOOST_CHECK((m_D_slice_1_u.at<0, 0>().value() == 1.0));
    BOOST_CHECK((m_D_slice_1_u.at<0, 1>().value() == 2.0));
    BOOST_CHECK((m_D_slice_1_u.at<1, 0>().value() == 5.0));
    BOOST_CHECK((m_D_slice_1_u.at<1, 1>().value() == 6.0));

    m_D_slice_1_u = bub::slice<rows_type_1, columns_type_1>(m_D1_u_d);
    BOOST_CHECK((m_D_slice_1_u.at<0, 0>().value() == 1.0));
    BOOST_CHECK((m_D_slice_1_u.at<0, 1>().value() == 2.0));
    BOOST_CHECK((m_D_slice_1_u.at<1, 0>().value() == 5.0));
    BOOST_CHECK((m_D_slice_1_u.at<1, 1>().value() == 6.0));

    typedef bub::result_of::slice<
        D_matrix_4x4_units_type,
        rows_type_2,
        columns_type_2
    >::type D_units_slice_type_2;
    D_units_slice_type_2 m_D_slice_2_u;
    m_D_slice_2_u = bub::slice<rows_type_2, columns_type_2>(m_D1_u);
    BOOST_CHECK((m_D_slice_2_u.at<0, 0>().value() == 3.0));
    BOOST_CHECK((m_D_slice_2_u.at<1, 0>().value() == 15.0));
    BOOST_CHECK((m_D_slice_2_u.at<2, 0>().value() == 11.0));

    m_D_slice_2_u = bub::slice<rows_type_2, columns_type_2>(m_D1_u_d);
    BOOST_CHECK((m_D_slice_2_u.at<0, 0>().value() == 3.0));
    BOOST_CHECK((m_D_slice_2_u.at<1, 0>().value() == 15.0));
    BOOST_CHECK((m_D_slice_2_u.at<2, 0>().value() == 11.0));

    typedef bub::result_of::slice<
        D_matrix_4x4_units_type,
        rows_type_3,
        columns_type_3
    >::type D_units_slice_type_3;
    D_units_slice_type_3 m_D_slice_3_u;
    m_D_slice_3_u = bub::slice<rows_type_3, columns_type_3>(m_D1_u);
    BOOST_CHECK((m_D_slice_3_u.at<0, 0>().value() == 2.0));
    BOOST_CHECK((m_D_slice_3_u.at<0, 1>().value() == 4.0));

    m_D_slice_3_u = bub::slice<rows_type_3, columns_type_3>(m_D1_u_d);
    BOOST_CHECK((m_D_slice_3_u.at<0, 0>().value() == 2.0));
    BOOST_CHECK((m_D_slice_3_u.at<0, 1>().value() == 4.0));

    return 0;
}
