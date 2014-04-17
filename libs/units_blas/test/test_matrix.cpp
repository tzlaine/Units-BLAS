// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

#include <boost/test/minimal.hpp>


namespace bub = boost::units_blas;

typedef bub::matrix<
    std::tuple<float>,
    std::tuple<int>
> A_matrix_2x1_float_int_type;

typedef bub::matrix<
    std::tuple<int>,
    std::tuple<int>
> B_matrix_2x1_int_type;

typedef bub::matrix<
    std::tuple<double, double>
> C_matrix_1x2_double_type;

typedef boost::units::quantity<boost::units::si::time> time_;
typedef boost::units::quantity<boost::units::si::length> length;
typedef boost::units::quantity<boost::units::si::dimensionless> dimensionless;

typedef bub::matrix<
    std::tuple<time_>,
    std::tuple<length>
> A_matrix_2x1_time_length_type;

typedef bub::matrix<
    std::tuple<dimensionless, dimensionless>
> B_matrix_2x1_length_type;

int test_main (int, char *[])
{
    // fundamental types

    int const i = 3;
    float const f = 2.0e8;

    A_matrix_2x1_float_int_type m_A1;
    m_A1.at<0, 0>() = f;
    m_A1.at<1, 0>() = i;

    A_matrix_2x1_float_int_type m_A2(m_A1);
    BOOST_CHECK((m_A2.at<0, 0>() == f));
    BOOST_CHECK((m_A2.at<1, 0>() == i));

    A_matrix_2x1_float_int_type m_A3;
    m_A3 = m_A1;
    BOOST_CHECK((m_A3.at<0, 0>() == f));
    BOOST_CHECK((m_A3.at<1, 0>() == i));

    B_matrix_2x1_int_type m_B1(m_A1);
    BOOST_CHECK((m_B1.at<0, 0>() == static_cast<int>(f)));
    BOOST_CHECK((m_B1.at<1, 0>() == i));

    B_matrix_2x1_int_type m_B2;
    m_B2 = m_A1;
    BOOST_CHECK((m_B2.at<0, 0>() == static_cast<int>(f)));
    BOOST_CHECK((m_B2.at<1, 0>() == i));

    BOOST_CHECK((m_A1.size() == 2));
    BOOST_CHECK((m_B1.size() == 2));
    BOOST_CHECK((m_A1.rows() == 2));
    BOOST_CHECK((m_B1.rows() == 2));
    BOOST_CHECK((m_A1.columns() == 1));
    BOOST_CHECK((m_B1.columns() == 1));

    A_matrix_2x1_float_int_type const m_A1_c(m_A1);
    BOOST_CHECK((m_A1_c.at<0, 0>() == f));
    BOOST_CHECK((m_A1_c.at<1, 0>() == i));

    A_matrix_2x1_float_int_type m_A4;
    m_A4.at<0, 0>() = f;
    m_A4.at<1, 0>() = i;
    B_matrix_2x1_int_type m_B3(m_A4);

    m_B3 += m_A4;
    BOOST_CHECK((m_B3.at<0, 0>() == static_cast<int>(static_cast<int>(f) + f)));
    BOOST_CHECK((m_B3.at<1, 0>() == i + i));

    m_B3 -= m_A4;
    BOOST_CHECK((m_B3.at<0, 0>() == static_cast<int>(static_cast<int>(static_cast<int>(f) + f) - static_cast<int>(f))));
    BOOST_CHECK((m_B3.at<1, 0>() == i + i - i));

    m_B3 = m_A4;
    m_B3 *= 2.0;
    BOOST_CHECK((m_B3.at<0, 0>() == static_cast<int>(static_cast<int>(f) * 2.0)));
    BOOST_CHECK((m_B3.at<1, 0>() == static_cast<int>(i * 2.0)));

    m_B3 = m_A4;
    m_B3 /= 2;
    BOOST_CHECK((m_B3.at<0, 0>() == static_cast<int>(f) / 2));
    BOOST_CHECK((m_B3.at<1, 0>() == i / 2));


    // unit types

    time_ const t = time_::from_value(1.0);
    length const l = length::from_value(2.0e8);

    A_matrix_2x1_time_length_type m_A1_u;
    m_A1_u.at<0, 0>() = t;
    m_A1_u.at<1, 0>() = l;

    A_matrix_2x1_time_length_type m_A2_u(m_A1_u);
    BOOST_CHECK((m_A2_u.at<0, 0>() == t));
    BOOST_CHECK((m_A2_u.at<1, 0>() == l));

    A_matrix_2x1_time_length_type m_A3_u;
    m_A3_u = m_A1_u;
    BOOST_CHECK((m_A3_u.at<0, 0>() == t));
    BOOST_CHECK((m_A3_u.at<1, 0>() == l));

    double const d1 = 1.0;
    double const d2 = 2.0;
    C_matrix_1x2_double_type m_C1_u;
    m_C1_u.at<0, 0>() = d1;
    m_C1_u.at<0, 1>() = d2;
    B_matrix_2x1_length_type m_B1_u(m_C1_u);
    BOOST_CHECK((m_B1_u.at<0, 0>() == static_cast<dimensionless>(d1)));
    BOOST_CHECK((m_B1_u.at<0, 1>() == static_cast<dimensionless>(d2)));

    BOOST_CHECK((m_A1_u.size() == 2));
    BOOST_CHECK((m_B1_u.size() == 2));
    BOOST_CHECK((m_A1_u.rows() == 2));
    BOOST_CHECK((m_B1_u.rows() == 1));
    BOOST_CHECK((m_A1_u.columns() == 1));
    BOOST_CHECK((m_B1_u.columns() == 2));

    A_matrix_2x1_time_length_type const m_A1_u_c(m_A1_u);
    BOOST_CHECK((m_A1_u_c.at<0, 0>() == t));
    BOOST_CHECK((m_A1_u_c.at<1, 0>() == l));

    A_matrix_2x1_time_length_type m_A4_u;
    m_A4_u.at<0, 0>() = t;
    m_A4_u.at<1, 0>() = l;
    A_matrix_2x1_time_length_type m_A5_u(m_A4_u);

    m_A5_u += m_A4_u;
    BOOST_CHECK((m_A5_u.at<0, 0>() == t + t));
    BOOST_CHECK((m_A5_u.at<1, 0>() == l + l));

    m_A5_u -= m_A4_u;
    BOOST_CHECK((m_A5_u.at<0, 0>() == t + t - t));
    BOOST_CHECK((m_A5_u.at<1, 0>() == l + l - l));

    m_A5_u = m_A4_u;
    m_A5_u *= 2.0;
    BOOST_CHECK((m_A5_u.at<0, 0>() == t * 2.0));
    BOOST_CHECK((m_A5_u.at<1, 0>() == l * 2.0));

    m_A5_u = m_A4_u;
    m_A5_u /= 2.0;
    BOOST_CHECK((m_A5_u.at<0, 0>() == t / 2.0));
    BOOST_CHECK((m_A5_u.at<1, 0>() == l / 2.0));

    return 0;
}
