// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "operations_tests.hpp"

#include <boost/test/minimal.hpp>


namespace bub = boost::units_blas;

struct derived_from_A_matrix_3x1_fundamentals_type :
    A_matrix_3x1_fundamentals_type
{};

struct derived_from_B_matrix_1x3_fundamentals_type :
    B_matrix_1x3_fundamentals_type
{};

int test_main (int, char *[])
{
    // funadamental types

    A_matrix_3x1_fundamentals_type m_A1;
    m_A1.at<0, 0>() = 1.0;
    m_A1.at<1, 0>() = 2.0;
    m_A1.at<2, 0>() = 3.0;

    derived_from_A_matrix_3x1_fundamentals_type m_A1_d;
    m_A1_d.at<0, 0>() = 1.0;
    m_A1_d.at<1, 0>() = 2.0;
    m_A1_d.at<2, 0>() = 3.0;

    A_matrix_3x1_fundamentals_type m_A2;
    m_A2.at<0, 0>() = 3.0;
    m_A2.at<1, 0>() = 2.0;
    m_A2.at<2, 0>() = 1.0;

    derived_from_A_matrix_3x1_fundamentals_type m_A2_d;
    m_A2_d.at<0, 0>() = 3.0;
    m_A2_d.at<1, 0>() = 2.0;
    m_A2_d.at<2, 0>() = 1.0;

    swap(m_A1, m_A2);
    BOOST_CHECK((m_A1.at<0, 0>() == 3.0));
    BOOST_CHECK((m_A1.at<1, 0>() == 2.0));
    BOOST_CHECK((m_A1.at<2, 0>() == 1.0));
    BOOST_CHECK((m_A2.at<0, 0>() == 1.0));
    BOOST_CHECK((m_A2.at<1, 0>() == 2.0));
    BOOST_CHECK((m_A2.at<2, 0>() == 3.0));
    swap(m_A1, m_A2);

    swap(m_A1_d, m_A2);
    BOOST_CHECK((m_A1_d.at<0, 0>() == 3.0));
    BOOST_CHECK((m_A1_d.at<1, 0>() == 2.0));
    BOOST_CHECK((m_A1_d.at<2, 0>() == 1.0));
    BOOST_CHECK((m_A2.at<0, 0>() == 1.0));
    BOOST_CHECK((m_A2.at<1, 0>() == 2.0));
    BOOST_CHECK((m_A2.at<2, 0>() == 3.0));
    swap(m_A1_d, m_A2);

    swap(m_A1, m_A2_d);
    BOOST_CHECK((m_A1.at<0, 0>() == 3.0));
    BOOST_CHECK((m_A1.at<1, 0>() == 2.0));
    BOOST_CHECK((m_A1.at<2, 0>() == 1.0));
    BOOST_CHECK((m_A2_d.at<0, 0>() == 1.0));
    BOOST_CHECK((m_A2_d.at<1, 0>() == 2.0));
    BOOST_CHECK((m_A2_d.at<2, 0>() == 3.0));
    swap(m_A1, m_A2_d);

    swap(m_A1_d, m_A2_d);
    BOOST_CHECK((m_A1_d.at<0, 0>() == 3.0));
    BOOST_CHECK((m_A1_d.at<1, 0>() == 2.0));
    BOOST_CHECK((m_A1_d.at<2, 0>() == 1.0));
    BOOST_CHECK((m_A2_d.at<0, 0>() == 1.0));
    BOOST_CHECK((m_A2_d.at<1, 0>() == 2.0));
    BOOST_CHECK((m_A2_d.at<2, 0>() == 3.0));
    swap(m_A1_d, m_A2_d);

    A_matrix_3x1_fundamentals_type m_A1_plus_A2;
    m_A1_plus_A2 = sum(m_A1, m_A2);
    BOOST_CHECK((m_A1_plus_A2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_A1_plus_A2.at<1, 0>() == 2.0 + 2.0));
    BOOST_CHECK((m_A1_plus_A2.at<2, 0>() == 3.0 + 1.0));
    m_A1_plus_A2 = m_A1 + m_A2;
    BOOST_CHECK((m_A1_plus_A2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_A1_plus_A2.at<1, 0>() == 2.0 + 2.0));
    BOOST_CHECK((m_A1_plus_A2.at<2, 0>() == 3.0 + 1.0));

    m_A1_plus_A2 = sum(m_A1_d, m_A2);
    BOOST_CHECK((m_A1_plus_A2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_A1_plus_A2.at<1, 0>() == 2.0 + 2.0));
    BOOST_CHECK((m_A1_plus_A2.at<2, 0>() == 3.0 + 1.0));
    m_A1_plus_A2 = m_A1_d + m_A2;
    BOOST_CHECK((m_A1_plus_A2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_A1_plus_A2.at<1, 0>() == 2.0 + 2.0));
    BOOST_CHECK((m_A1_plus_A2.at<2, 0>() == 3.0 + 1.0));

    m_A1_plus_A2 = sum(m_A1, m_A2_d);
    BOOST_CHECK((m_A1_plus_A2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_A1_plus_A2.at<1, 0>() == 2.0 + 2.0));
    BOOST_CHECK((m_A1_plus_A2.at<2, 0>() == 3.0 + 1.0));
    m_A1_plus_A2 = m_A1 + m_A2_d;
    BOOST_CHECK((m_A1_plus_A2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_A1_plus_A2.at<1, 0>() == 2.0 + 2.0));
    BOOST_CHECK((m_A1_plus_A2.at<2, 0>() == 3.0 + 1.0));

    m_A1_plus_A2 = sum(m_A1_d, m_A2_d);
    BOOST_CHECK((m_A1_plus_A2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_A1_plus_A2.at<1, 0>() == 2.0 + 2.0));
    BOOST_CHECK((m_A1_plus_A2.at<2, 0>() == 3.0 + 1.0));
    m_A1_plus_A2 = m_A1_d + m_A2_d;
    BOOST_CHECK((m_A1_plus_A2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_A1_plus_A2.at<1, 0>() == 2.0 + 2.0));
    BOOST_CHECK((m_A1_plus_A2.at<2, 0>() == 3.0 + 1.0));

    A_matrix_3x1_fundamentals_type m_A1_minus_A2;
    m_A1_minus_A2 = diff(m_A1, m_A2);
    BOOST_CHECK((m_A1_minus_A2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_A1_minus_A2.at<1, 0>() == 2.0 - 2.0));
    BOOST_CHECK((m_A1_minus_A2.at<2, 0>() == 3.0 - 1.0));
    m_A1_minus_A2 = m_A1 - m_A2;
    BOOST_CHECK((m_A1_minus_A2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_A1_minus_A2.at<1, 0>() == 2.0 - 2.0));
    BOOST_CHECK((m_A1_minus_A2.at<2, 0>() == 3.0 - 1.0));

    m_A1_minus_A2 = diff(m_A1_d, m_A2);
    BOOST_CHECK((m_A1_minus_A2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_A1_minus_A2.at<1, 0>() == 2.0 - 2.0));
    BOOST_CHECK((m_A1_minus_A2.at<2, 0>() == 3.0 - 1.0));
    m_A1_minus_A2 = m_A1_d - m_A2;
    BOOST_CHECK((m_A1_minus_A2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_A1_minus_A2.at<1, 0>() == 2.0 - 2.0));
    BOOST_CHECK((m_A1_minus_A2.at<2, 0>() == 3.0 - 1.0));

    m_A1_minus_A2 = diff(m_A1, m_A2_d);
    BOOST_CHECK((m_A1_minus_A2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_A1_minus_A2.at<1, 0>() == 2.0 - 2.0));
    BOOST_CHECK((m_A1_minus_A2.at<2, 0>() == 3.0 - 1.0));
    m_A1_minus_A2 = m_A1 - m_A2_d;
    BOOST_CHECK((m_A1_minus_A2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_A1_minus_A2.at<1, 0>() == 2.0 - 2.0));
    BOOST_CHECK((m_A1_minus_A2.at<2, 0>() == 3.0 - 1.0));

    m_A1_minus_A2 = diff(m_A1_d, m_A2_d);
    BOOST_CHECK((m_A1_minus_A2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_A1_minus_A2.at<1, 0>() == 2.0 - 2.0));
    BOOST_CHECK((m_A1_minus_A2.at<2, 0>() == 3.0 - 1.0));
    m_A1_minus_A2 = m_A1_d - m_A2_d;
    BOOST_CHECK((m_A1_minus_A2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_A1_minus_A2.at<1, 0>() == 2.0 - 2.0));
    BOOST_CHECK((m_A1_minus_A2.at<2, 0>() == 3.0 - 1.0));

    B_matrix_1x3_fundamentals_type m_B1;
    m_B1.at<0, 0>() = 1.0;
    m_B1.at<0, 1>() = 2.0;
    m_B1.at<0, 2>() = 3.0;

    derived_from_B_matrix_1x3_fundamentals_type m_B1_d;
    m_B1_d.at<0, 0>() = 1.0;
    m_B1_d.at<0, 1>() = 2.0;
    m_B1_d.at<0, 2>() = 3.0;

    B_matrix_1x3_fundamentals_type m_B2;
    m_B2.at<0, 0>() = 3.0;
    m_B2.at<0, 1>() = 2.0;
    m_B2.at<0, 2>() = 1.0;

    derived_from_B_matrix_1x3_fundamentals_type m_B2_d;
    m_B2_d.at<0, 0>() = 3.0;
    m_B2_d.at<0, 1>() = 2.0;
    m_B2_d.at<0, 2>() = 1.0;

    swap(m_B1, m_B2);
    BOOST_CHECK((m_B1.at<0, 0>() == 3.0));
    BOOST_CHECK((m_B1.at<0, 1>() == 2.0));
    BOOST_CHECK((m_B1.at<0, 2>() == 1.0));
    BOOST_CHECK((m_B2.at<0, 0>() == 1.0));
    BOOST_CHECK((m_B2.at<0, 1>() == 2.0));
    BOOST_CHECK((m_B2.at<0, 2>() == 3.0));
    swap(m_B1, m_B2);

    B_matrix_1x3_fundamentals_type m_B1_plus_B2 = sum(m_B1, m_B2);
    BOOST_CHECK((m_B1_plus_B2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_B1_plus_B2.at<0, 1>() == 2.0 + 2.0));
    BOOST_CHECK((m_B1_plus_B2.at<0, 2>() == 3.0 + 1.0));
    m_B1_plus_B2 = m_B1 + m_B2;
    BOOST_CHECK((m_B1_plus_B2.at<0, 0>() == 1.0 + 3.0));
    BOOST_CHECK((m_B1_plus_B2.at<0, 1>() == 2.0 + 2.0));
    BOOST_CHECK((m_B1_plus_B2.at<0, 2>() == 3.0 + 1.0));

    B_matrix_1x3_fundamentals_type m_B1_minus_B2 = diff(m_B1, m_B2);
    BOOST_CHECK((m_B1_minus_B2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_B1_minus_B2.at<0, 1>() == 2.0 - 2.0));
    BOOST_CHECK((m_B1_minus_B2.at<0, 2>() == 3.0 - 1.0));
    m_B1_minus_B2 = m_B1 - m_B2;
    BOOST_CHECK((m_B1_minus_B2.at<0, 0>() == 1.0 - 3.0));
    BOOST_CHECK((m_B1_minus_B2.at<0, 1>() == 2.0 - 2.0));
    BOOST_CHECK((m_B1_minus_B2.at<0, 2>() == 3.0 - 1.0));

    typedef bub::result_of::matrix_product<
        A_matrix_3x1_fundamentals_type,
        B_matrix_1x3_fundamentals_type
    >::type AxB_type;
    AxB_type m_A1_times_B1;
    m_A1_times_B1 = prod(m_A1, m_B1);
    BOOST_CHECK((m_A1_times_B1.size() == 9));
    BOOST_CHECK((m_A1_times_B1.at<0, 0>() == 1.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 1>() == 1.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 2>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 0>() == 2.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 2>() == 2.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 0>() == 3.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 1>() == 3.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 2>() == 3.0 * 3.0));
    m_A1_times_B1 = m_A1 * m_B1;
    BOOST_CHECK((m_A1_times_B1.at<0, 0>() == 1.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 1>() == 1.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 2>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 0>() == 2.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 2>() == 2.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 0>() == 3.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 1>() == 3.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 2>() == 3.0 * 3.0));

    m_A1_times_B1 = prod(m_A1_d, m_B1);
    BOOST_CHECK((m_A1_times_B1.size() == 9));
    BOOST_CHECK((m_A1_times_B1.at<0, 0>() == 1.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 1>() == 1.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 2>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 0>() == 2.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 2>() == 2.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 0>() == 3.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 1>() == 3.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 2>() == 3.0 * 3.0));
    m_A1_times_B1 = m_A1_d * m_B1;
    BOOST_CHECK((m_A1_times_B1.at<0, 0>() == 1.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 1>() == 1.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 2>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 0>() == 2.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 2>() == 2.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 0>() == 3.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 1>() == 3.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 2>() == 3.0 * 3.0));

    m_A1_times_B1 = prod(m_A1, m_B1_d);
    BOOST_CHECK((m_A1_times_B1.size() == 9));
    BOOST_CHECK((m_A1_times_B1.at<0, 0>() == 1.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 1>() == 1.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 2>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 0>() == 2.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 2>() == 2.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 0>() == 3.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 1>() == 3.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 2>() == 3.0 * 3.0));
    m_A1_times_B1 = m_A1 * m_B1_d;
    BOOST_CHECK((m_A1_times_B1.at<0, 0>() == 1.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 1>() == 1.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 2>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 0>() == 2.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 2>() == 2.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 0>() == 3.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 1>() == 3.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 2>() == 3.0 * 3.0));

    m_A1_times_B1 = prod(m_A1_d, m_B1_d);
    BOOST_CHECK((m_A1_times_B1.size() == 9));
    BOOST_CHECK((m_A1_times_B1.at<0, 0>() == 1.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 1>() == 1.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 2>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 0>() == 2.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 2>() == 2.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 0>() == 3.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 1>() == 3.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 2>() == 3.0 * 3.0));
    m_A1_times_B1 = m_A1_d * m_B1_d;
    BOOST_CHECK((m_A1_times_B1.at<0, 0>() == 1.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 1>() == 1.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<0, 2>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 0>() == 2.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<1, 2>() == 2.0 * 3.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 0>() == 3.0 * 1.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 1>() == 3.0 * 2.0));
    BOOST_CHECK((m_A1_times_B1.at<2, 2>() == 3.0 * 3.0));

    typedef bub::result_of::matrix_product<
        B_matrix_1x3_fundamentals_type,
        A_matrix_3x1_fundamentals_type
    >::type BxA_type;
    BxA_type m_B1_times_A1;
    m_B1_times_A1 = prod(m_B1, m_A1);
    BOOST_CHECK((m_B1_times_A1.size() == 1));
    BOOST_CHECK((m_B1_times_A1.at<0, 0>() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));
    m_B1_times_A1 = m_B1 * m_A1;
    BOOST_CHECK((m_B1_times_A1.at<0, 0>() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));

    m_B1_times_A1 = prod(m_B1_d, m_A1);
    BOOST_CHECK((m_B1_times_A1.size() == 1));
    BOOST_CHECK((m_B1_times_A1.at<0, 0>() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));
    m_B1_times_A1 = m_B1_d * m_A1;
    BOOST_CHECK((m_B1_times_A1.at<0, 0>() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));

    m_B1_times_A1 = prod(m_B1, m_A1_d);
    BOOST_CHECK((m_B1_times_A1.size() == 1));
    BOOST_CHECK((m_B1_times_A1.at<0, 0>() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));
    m_B1_times_A1 = m_B1 * m_A1_d;
    BOOST_CHECK((m_B1_times_A1.at<0, 0>() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));

    m_B1_times_A1 = prod(m_B1_d, m_A1_d);
    BOOST_CHECK((m_B1_times_A1.size() == 1));
    BOOST_CHECK((m_B1_times_A1.at<0, 0>() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));
    m_B1_times_A1 = m_B1_d * m_A1_d;
    BOOST_CHECK((m_B1_times_A1.at<0, 0>() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));

    E_matrix_1x1_fundamentals_type m_E;
    m_E.at<0, 0>() = -2.0;
    typedef bub::result_of::matrix_product<
        E_matrix_1x1_fundamentals_type,
        E_matrix_1x1_fundamentals_type
    >::type ExE_fundamentals_type;
    ExE_fundamentals_type m_E_times_E = prod(m_E, m_E);
    BOOST_CHECK((m_E_times_E.at<0, 0>() == 4.0));
    m_E_times_E = m_E * m_E;
    BOOST_CHECK((m_E_times_E.at<0, 0>() == 4.0));

    A_matrix_3x1_fundamentals_type m_A1_elem_times_A2;
    m_A1_elem_times_A2 = element_prod(m_A1, m_A2);
    BOOST_CHECK((m_A1_elem_times_A2.at<0, 0>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_elem_times_A2.at<1, 0>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_elem_times_A2.at<2, 0>() == 3.0 * 1.0));
    A_matrix_3x1_fundamentals_type m_A1_elem_div_A2;
    m_A1_elem_div_A2 = element_div(m_A1, m_A2);
    BOOST_CHECK((m_A1_elem_div_A2.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_A1_elem_div_A2.at<1, 0>() == 2.0 / 2.0));
    BOOST_CHECK((m_A1_elem_div_A2.at<2, 0>() == 3.0 / 1.0));

    m_A1_elem_times_A2 = element_prod(m_A1_d, m_A2);
    BOOST_CHECK((m_A1_elem_times_A2.at<0, 0>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_elem_times_A2.at<1, 0>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_elem_times_A2.at<2, 0>() == 3.0 * 1.0));
    m_A1_elem_div_A2 = element_div(m_A1_d, m_A2);
    BOOST_CHECK((m_A1_elem_div_A2.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_A1_elem_div_A2.at<1, 0>() == 2.0 / 2.0));
    BOOST_CHECK((m_A1_elem_div_A2.at<2, 0>() == 3.0 / 1.0));

    m_A1_elem_times_A2 = element_prod(m_A1, m_A2_d);
    BOOST_CHECK((m_A1_elem_times_A2.at<0, 0>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_elem_times_A2.at<1, 0>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_elem_times_A2.at<2, 0>() == 3.0 * 1.0));
    m_A1_elem_div_A2 = element_div(m_A1, m_A2_d);
    BOOST_CHECK((m_A1_elem_div_A2.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_A1_elem_div_A2.at<1, 0>() == 2.0 / 2.0));
    BOOST_CHECK((m_A1_elem_div_A2.at<2, 0>() == 3.0 / 1.0));

    m_A1_elem_times_A2 = element_prod(m_A1_d, m_A2_d);
    BOOST_CHECK((m_A1_elem_times_A2.at<0, 0>() == 1.0 * 3.0));
    BOOST_CHECK((m_A1_elem_times_A2.at<1, 0>() == 2.0 * 2.0));
    BOOST_CHECK((m_A1_elem_times_A2.at<2, 0>() == 3.0 * 1.0));
    m_A1_elem_div_A2 = element_div(m_A1_d, m_A2_d);
    BOOST_CHECK((m_A1_elem_div_A2.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_A1_elem_div_A2.at<1, 0>() == 2.0 / 2.0));
    BOOST_CHECK((m_A1_elem_div_A2.at<2, 0>() == 3.0 / 1.0));

    B_matrix_1x3_fundamentals_type m_B1_elem_times_B2;
    m_B1_elem_times_B2 = element_prod(m_B1, m_B2);
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 0>() == 1.0 * 3.0));
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 2>() == 3.0 * 1.0));
    B_matrix_1x3_fundamentals_type m_B1_elem_div_B2;
    m_B1_elem_div_B2 = element_div(m_B1, m_B2);
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 1>() == 2.0 / 2.0));
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 2>() == 3.0 / 1.0));

    m_B1_elem_times_B2 = element_prod(m_B1_d, m_B2);
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 0>() == 1.0 * 3.0));
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 2>() == 3.0 * 1.0));
    m_B1_elem_div_B2 = element_div(m_B1_d, m_B2);
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 1>() == 2.0 / 2.0));
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 2>() == 3.0 / 1.0));

    m_B1_elem_times_B2 = element_prod(m_B1, m_B2_d);
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 0>() == 1.0 * 3.0));
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 2>() == 3.0 * 1.0));
    m_B1_elem_div_B2 = element_div(m_B1, m_B2_d);
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 1>() == 2.0 / 2.0));
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 2>() == 3.0 / 1.0));

    m_B1_elem_times_B2 = element_prod(m_B1_d, m_B2_d);
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 0>() == 1.0 * 3.0));
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 1>() == 2.0 * 2.0));
    BOOST_CHECK((m_B1_elem_times_B2.at<0, 2>() == 3.0 * 1.0));
    m_B1_elem_div_B2 = element_div(m_B1_d, m_B2_d);
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 0>() == 1.0 / 3.0));
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 1>() == 2.0 / 2.0));
    BOOST_CHECK((m_B1_elem_div_B2.at<0, 2>() == 3.0 / 1.0));


    // unit types

    C_matrix_3x1_units_type m_C1_u;
    m_C1_u.at<0, 0>() = dimensionless::from_value(1.0);
    m_C1_u.at<1, 0>() = dimensionless::from_value(2.0);
    m_C1_u.at<2, 0>() = dimensionless::from_value(3.0);

    C_matrix_3x1_units_type m_C2_u;
    m_C2_u.at<0, 0>() = dimensionless::from_value(3.0);
    m_C2_u.at<1, 0>() = dimensionless::from_value(2.0);
    m_C2_u.at<2, 0>() = dimensionless::from_value(1.0);

    swap(m_C1_u, m_C2_u);
    BOOST_CHECK((m_C1_u.at<0, 0>().value() == 3.0));
    BOOST_CHECK((m_C1_u.at<1, 0>().value() == 2.0));
    BOOST_CHECK((m_C1_u.at<2, 0>().value() == 1.0));
    BOOST_CHECK((m_C2_u.at<0, 0>().value() == 1.0));
    BOOST_CHECK((m_C2_u.at<1, 0>().value() == 2.0));
    BOOST_CHECK((m_C2_u.at<2, 0>().value() == 3.0));
    swap(m_C1_u, m_C2_u);

    C_matrix_3x1_units_type m_C1_u_plus_C2_u = sum(m_C1_u, m_C2_u);
    BOOST_CHECK((m_C1_u_plus_C2_u.at<0, 0>().value() == 1.0 + 3.0));
    BOOST_CHECK((m_C1_u_plus_C2_u.at<1, 0>().value() == 2.0 + 2.0));
    BOOST_CHECK((m_C1_u_plus_C2_u.at<2, 0>().value() == 3.0 + 1.0));
    m_C1_u_plus_C2_u = m_C1_u + m_C2_u;
    BOOST_CHECK((m_C1_u_plus_C2_u.at<0, 0>().value() == 1.0 + 3.0));
    BOOST_CHECK((m_C1_u_plus_C2_u.at<1, 0>().value() == 2.0 + 2.0));
    BOOST_CHECK((m_C1_u_plus_C2_u.at<2, 0>().value() == 3.0 + 1.0));

    C_matrix_3x1_units_type m_C1_u_minus_C2_u = diff(m_C1_u, m_C2_u);
    BOOST_CHECK((m_C1_u_minus_C2_u.at<0, 0>().value() == 1.0 - 3.0));
    BOOST_CHECK((m_C1_u_minus_C2_u.at<1, 0>().value() == 2.0 - 2.0));
    BOOST_CHECK((m_C1_u_minus_C2_u.at<2, 0>().value() == 3.0 - 1.0));
    m_C1_u_minus_C2_u = m_C1_u - m_C2_u;
    BOOST_CHECK((m_C1_u_minus_C2_u.at<0, 0>().value() == 1.0 - 3.0));
    BOOST_CHECK((m_C1_u_minus_C2_u.at<1, 0>().value() == 2.0 - 2.0));
    BOOST_CHECK((m_C1_u_minus_C2_u.at<2, 0>().value() == 3.0 - 1.0));

    B_matrix_1x3_units_type m_B1_u;
    m_B1_u.at<0, 0>() = time_::from_value(1.0);
    m_B1_u.at<0, 1>() = time_::from_value(2.0);
    m_B1_u.at<0, 2>() = time_::from_value(3.0);

    B_matrix_1x3_units_type m_B2_u;
    m_B2_u.at<0, 0>() = time_::from_value(3.0);
    m_B2_u.at<0, 1>() = time_::from_value(2.0);
    m_B2_u.at<0, 2>() = time_::from_value(1.0);

    swap(m_B1_u, m_B2_u);
    BOOST_CHECK((m_B1_u.at<0, 0>().value() == 3.0));
    BOOST_CHECK((m_B1_u.at<0, 1>().value() == 2.0));
    BOOST_CHECK((m_B1_u.at<0, 2>().value() == 1.0));
    BOOST_CHECK((m_B2_u.at<0, 0>().value() == 1.0));
    BOOST_CHECK((m_B2_u.at<0, 1>().value() == 2.0));
    BOOST_CHECK((m_B2_u.at<0, 2>().value() == 3.0));
    swap(m_B1_u, m_B2_u);

    B_matrix_1x3_units_type m_B1_u_plus_B2_u = sum(m_B1_u, m_B2_u);
    BOOST_CHECK((m_B1_u_plus_B2_u.at<0, 0>().value() == 1.0 + 3.0));
    BOOST_CHECK((m_B1_u_plus_B2_u.at<0, 1>().value() == 2.0 + 2.0));
    BOOST_CHECK((m_B1_u_plus_B2_u.at<0, 2>().value() == 3.0 + 1.0));
    m_B1_u_plus_B2_u = m_B1_u + m_B2_u;
    BOOST_CHECK((m_B1_u_plus_B2_u.at<0, 0>().value() == 1.0 + 3.0));
    BOOST_CHECK((m_B1_u_plus_B2_u.at<0, 1>().value() == 2.0 + 2.0));
    BOOST_CHECK((m_B1_u_plus_B2_u.at<0, 2>().value() == 3.0 + 1.0));

    B_matrix_1x3_units_type m_B1_u_minus_B2_u = diff(m_B1_u, m_B2_u);
    BOOST_CHECK((m_B1_u_minus_B2_u.at<0, 0>().value() == 1.0 - 3.0));
    BOOST_CHECK((m_B1_u_minus_B2_u.at<0, 1>().value() == 2.0 - 2.0));
    BOOST_CHECK((m_B1_u_minus_B2_u.at<0, 2>().value() == 3.0 - 1.0));
    m_B1_u_minus_B2_u = m_B1_u - m_B2_u;
    BOOST_CHECK((m_B1_u_minus_B2_u.at<0, 0>().value() == 1.0 - 3.0));
    BOOST_CHECK((m_B1_u_minus_B2_u.at<0, 1>().value() == 2.0 - 2.0));
    BOOST_CHECK((m_B1_u_minus_B2_u.at<0, 2>().value() == 3.0 - 1.0));

    typedef bub::result_of::matrix_product<
        C_matrix_3x1_units_type,
        B_matrix_1x3_units_type
    >::type CxB_type;
    CxB_type m_C1_u_times_B1_u = prod(m_C1_u, m_B1_u);
    BOOST_CHECK((m_C1_u_times_B1_u.size() == 9));
    BOOST_CHECK((m_C1_u_times_B1_u.at<0, 0>().value() == 1.0 * 1.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<0, 1>().value() == 1.0 * 2.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<0, 2>().value() == 1.0 * 3.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<1, 0>().value() == 2.0 * 1.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<1, 1>().value() == 2.0 * 2.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<1, 2>().value() == 2.0 * 3.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<2, 0>().value() == 3.0 * 1.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<2, 1>().value() == 3.0 * 2.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<2, 2>().value() == 3.0 * 3.0));
    m_C1_u_times_B1_u = m_C1_u * m_B1_u;
    BOOST_CHECK((m_C1_u_times_B1_u.at<0, 0>().value() == 1.0 * 1.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<0, 1>().value() == 1.0 * 2.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<0, 2>().value() == 1.0 * 3.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<1, 0>().value() == 2.0 * 1.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<1, 1>().value() == 2.0 * 2.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<1, 2>().value() == 2.0 * 3.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<2, 0>().value() == 3.0 * 1.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<2, 1>().value() == 3.0 * 2.0));
    BOOST_CHECK((m_C1_u_times_B1_u.at<2, 2>().value() == 3.0 * 3.0));

    E_matrix_1x1_units_type m_E_u;
    m_E_u.at<0, 0>() = length::from_value(-2.0);
    typedef bub::result_of::matrix_product<
        E_matrix_1x1_units_type,
        E_matrix_1x1_units_type
    >::type ExE_units_type;
    ExE_units_type m_E_times_E_u = prod(m_E_u, m_E_u);
    BOOST_CHECK((m_E_times_E_u.at<0, 0>().value() == 4.0));
    m_E_times_E_u = m_E_u * m_E_u;
    BOOST_CHECK((m_E_times_E_u.at<0, 0>().value() == 4.0));

    typedef bub::result_of::matrix_product<
        B_matrix_1x3_units_type,
        C_matrix_3x1_units_type
    >::type BxC_type;
    BxC_type m_B1_u_times_C1_u = prod(m_B1_u, m_C1_u);
    BOOST_CHECK((m_B1_u_times_C1_u.size() == 1));
    BOOST_CHECK((m_B1_u_times_C1_u.at<0, 0>().value() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));
    m_B1_u_times_C1_u = m_B1_u * m_C1_u;
    BOOST_CHECK((m_B1_u_times_C1_u.at<0, 0>().value() == 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0));

    typedef bub::result_of::matrix_element_product<
        C_matrix_3x1_units_type,
        C_matrix_3x1_units_type
    >::type CxC_type;
    CxC_type m_C1_u_elem_times_C2_u = element_prod(m_C1_u, m_C2_u);
    BOOST_CHECK((m_C1_u_elem_times_C2_u.at<0, 0>().value() == 1.0 * 3.0));
    BOOST_CHECK((m_C1_u_elem_times_C2_u.at<1, 0>().value() == 2.0 * 2.0));
    BOOST_CHECK((m_C1_u_elem_times_C2_u.at<2, 0>().value() == 3.0 * 1.0));
    typedef bub::result_of::matrix_element_quotient<
        C_matrix_3x1_units_type,
        C_matrix_3x1_units_type
    >::type CdC_type;
    CdC_type m_C1_u_elem_div_C2_u = element_div(m_C1_u, m_C2_u);
    BOOST_CHECK((m_C1_u_elem_div_C2_u.at<0, 0>().value() == 1.0 / 3.0));
    BOOST_CHECK((m_C1_u_elem_div_C2_u.at<1, 0>().value() == 2.0 / 2.0));
    BOOST_CHECK((m_C1_u_elem_div_C2_u.at<2, 0>().value() == 3.0 / 1.0));

    typedef bub::result_of::matrix_element_product<
        B_matrix_1x3_units_type,
        B_matrix_1x3_units_type
    >::type BxB_type;
    BxB_type m_B1_u_elem_times_B2_u = element_prod(m_B1_u, m_B2_u);
    BOOST_CHECK((m_B1_u_elem_times_B2_u.at<0, 0>().value() == 1.0 * 3.0));
    BOOST_CHECK((m_B1_u_elem_times_B2_u.at<0, 1>().value() == 2.0 * 2.0));
    BOOST_CHECK((m_B1_u_elem_times_B2_u.at<0, 2>().value() == 3.0 * 1.0));
    typedef bub::result_of::matrix_element_quotient<
        B_matrix_1x3_units_type,
        B_matrix_1x3_units_type
    >::type BdB_type;
    BdB_type m_B1_u_elem_div_B2_u = element_div(m_B1_u, m_B2_u);
    BOOST_CHECK((m_B1_u_elem_div_B2_u.at<0, 0>().value() == 1.0 / 3.0));
    BOOST_CHECK((m_B1_u_elem_div_B2_u.at<0, 1>().value() == 2.0 / 2.0));
    BOOST_CHECK((m_B1_u_elem_div_B2_u.at<0, 2>().value() == 3.0 / 1.0));

    return 0;
}
