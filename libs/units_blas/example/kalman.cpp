//[units_blas_boilerplate
#include <boost/units_blas.hpp>

using namespace boost;
//]

//[boost_units_typedefs
#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/frequency.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

typedef units::quantity<units::SI::time> s;
typedef units::quantity<units::SI::length> m;
typedef units::quantity<units::SI::velocity> m_per_s;
typedef units::quantity<units::SI::dimensionless> none;
typedef units::quantity<units::SI::frequency> hz;

namespace boost { namespace units { namespace SI {
    typedef derived_dimension<length_base_dimension, 2>::type length_squared_dimension;
    typedef derived_dimension<length_base_dimension, 2, time_base_dimension, -1>::type length_squared_per_time_dimension;
    typedef derived_dimension<length_base_dimension, 2, time_base_dimension, -2>::type length_squared_per_time_squared_dimension;
    typedef unit<length_squared_dimension, system> length_squared;
    typedef unit<length_squared_per_time_dimension, system> length_squared_per_time;
    typedef unit<length_squared_per_time_squared_dimension, system> length_squared_per_time_squared;
} } }

typedef units::quantity<units::SI::length_squared> m2;
typedef units::quantity<units::SI::length_squared_per_time> m2_per_s;
typedef units::quantity<units::SI::length_squared_per_time_squared> m2_per_s2;
//]

//[meas_type
typedef units_blas::matrix<
    fusion::vector<
        fusion::vector<m>,
        fusion::vector<m>
    >
> meas_type;
//]

//[state_type
typedef units_blas::matrix<
    fusion::vector<
        fusion::vector<m>,
        fusion::vector<m_per_s>,
        fusion::vector<m>,
        fusion::vector<m_per_s>
    >
> state_type;
//]

//[P_type
typedef units_blas::matrix<
    fusion::vector<
        fusion::vector<m2,       m2_per_s,  m2,       m2_per_s>,
        fusion::vector<m2_per_s, m2_per_s2, m2_per_s, m2_per_s2>,
        fusion::vector<m2,       m2_per_s,  m2,       m2_per_s>,
        fusion::vector<m2_per_s, m2_per_s2, m2_per_s, m2_per_s2>
    >
> P_type;
//]

//[Q_type
typedef P_type Q_type;
//]

//[F_type
typedef units_blas::matrix<
    fusion::vector<
        fusion::vector<none, s,    none, s>,
        fusion::vector<hz,   none, hz,   none>,
        fusion::vector<none, s,    none, s>,
        fusion::vector<hz,   none, hz,   none>
    >
> F_type;
//]

//[R_type
typedef units_blas::matrix<
    fusion::vector<
        fusion::vector<m2, m2>,
        fusion::vector<m2, m2>
    >
> R_type;
//]

//[H_type
typedef units_blas::matrix<
    fusion::vector<
        fusion::vector<none, s, none, s>,
        fusion::vector<none, s, none, s>
    >
> H_type;
//]

//[S_type
typedef R_type S_type;
//]

//[K_type
typedef units_blas::matrix<
    fusion::vector<
        fusion::vector<none, none>,
        fusion::vector<hz,   hz>,
        fusion::vector<none, none>,
        fusion::vector<hz,   hz>
    >
> K_type;
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
