//[units_blas_boilerplate
#include <boost/units_blas.hpp>

#include <iostream>
#include <boost/units/io.hpp>

using namespace boost;
//]

//[boost_units_typedefs
#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/frequency.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

using s = units::quantity<units::si::time>;
using m = units::quantity<units::si::length>;
using m_per_s = units::quantity<units::si::velocity>;
using none_ = units::quantity<units::si::dimensionless>;
using hz = units::quantity<units::si::frequency>;

namespace boost { namespace units { namespace si {

    using length_squared_dimension =
        derived_dimension<length_base_dimension, 2>::type;
    using length_squared_per_time_dimension =
        derived_dimension<length_base_dimension, 2, time_base_dimension, -1>::type;
    using length_squared_per_time_squared_dimension =
        derived_dimension<length_base_dimension, 2, time_base_dimension, -2>::type;

    using length_squared =
        unit<length_squared_dimension, system>;
    using length_squared_per_time =
        unit<length_squared_per_time_dimension, system>;
    using length_squared_per_time_squared =
        unit<length_squared_per_time_squared_dimension, system>;

} } }

using m2 = units::quantity<units::si::length_squared>;
using m2_per_s = units::quantity<units::si::length_squared_per_time>;
using m2_per_s2 = units::quantity<units::si::length_squared_per_time_squared>;
//]

//[meas_type
using meas_type = units_blas::matrix<
    std::tuple<m>,
    std::tuple<m>
>;
//]

//[state_type
using state_type = units_blas::matrix<
    std::tuple<m>,
    std::tuple<m_per_s>,
    std::tuple<m>,
    std::tuple<m_per_s>
>;
//]

//[P_type
using P_type = units_blas::matrix<
    std::tuple<m2,       m2_per_s,  m2,       m2_per_s>,
    std::tuple<m2_per_s, m2_per_s2, m2_per_s, m2_per_s2>,
    std::tuple<m2,       m2_per_s,  m2,       m2_per_s>,
    std::tuple<m2_per_s, m2_per_s2, m2_per_s, m2_per_s2>
>;
//]

//[Q_type
using Q_type = P_type;
//]

//[F_type
using F_type = units_blas::matrix<
    std::tuple<none_, s,     none_, s>,
    std::tuple<hz,    none_, hz,    none_>,
    std::tuple<none_, s,     none_, s>,
    std::tuple<hz,    none_, hz,    none_>
>;
//]

//[R_type
using R_type = units_blas::matrix<
    std::tuple<m2, m2>,
    std::tuple<m2, m2>
>;
//]

//[H_type
using H_type = units_blas::matrix<
    std::tuple<none_, s, none_, s>,
    std::tuple<none_, s, none_, s>
>;
//]

//[S_type
using S_type = R_type;
//]

//[K_type
using K_type = units_blas::matrix<
    std::tuple<none_, none_>,
    std::tuple<hz,    hz>,
    std::tuple<none_, none_>,
    std::tuple<hz,    hz>
>;
//]

int main()
{
//[prediction_step_1
    // prediction step variables
    state_type x;
    P_type P;

    F_type F;
    Q_type Q;
//]

    x.at<0, 0>() = m::from_value(0.0);
    x.at<1, 0>() = m_per_s::from_value(1.0);
    x.at<2, 0>() = m::from_value(0.0);
    x.at<3, 0>() = m_per_s::from_value(1.0);

    P.at<0, 0>() = m2::from_value(1.0);
    P.at<1, 1>() = m2_per_s2::from_value(1.0);
    P.at<2, 2>() = m2::from_value(1.0);
    P.at<3, 3>() = m2_per_s2::from_value(1.0);

    F.at<0, 0>() = 1.0;
    F.at<1, 1>() = 1.0;
    F.at<2, 2>() = 1.0;
    F.at<3, 3>() = 1.0;

    Q = P;

//[prediction_step_2
    // prediction step
    state_type x_predict = F * x;
    P_type P_predict = F * P * transpose(F) + Q;
//]

//[update_step_1
    // update step variables
    H_type H;
    R_type R;
    meas_type z;
//]

    H.at<0, 0>() = 1.0;
    H.at<1, 2>() = 1.0;

    R.at<0, 0>() = m2::from_value(1.0);
    R.at<1, 1>() = m2::from_value(1.0);

    z.at<0, 0>() = m::from_value(1.0);
    z.at<1, 0>() = m::from_value(1.0);

//[update_step_2
    // update step
    using namespace boost::units_blas;
    S_type S = H * P_predict * transpose(H) + R;
    K_type K = P_predict * transpose(H) * inverse(S);
    x = x_predict + K * (z - H * x_predict);
    P = P_predict + K * S * transpose(K);
//]

    return 0;
}
