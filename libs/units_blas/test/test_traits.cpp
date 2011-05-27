// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/make_matrix.hpp>
#include <boost/units_blas/traits.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

#include <boost/test/minimal.hpp>


namespace bub = boost::units_blas;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<int>,
        boost::fusion::vector<long>,
        boost::fusion::vector<float>,
        boost::fusion::vector<double>
    >
>::type matrix_4x1_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, long, float, double>
    >
>::type matrix_1x4_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, long, float, double>,
        boost::fusion::vector<long, int, float, double>,
        boost::fusion::vector<int, float, long, double>,
        boost::fusion::vector<int, long, double, float>
    >
>::type matrix_4x4_type;

typedef boost::units::quantity<boost::units::si::time> time_;
typedef boost::units::quantity<boost::units::si::length> length;
typedef boost::units::quantity<boost::units::si::velocity> velocity;
typedef boost::units::quantity<boost::units::si::dimensionless> dimensionless;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_>,
        boost::fusion::vector<length>,
        boost::fusion::vector<velocity>,
        boost::fusion::vector<dimensionless>
    >
>::type matrix_4x1_units_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_, length, velocity, dimensionless>
    >
>::type matrix_1x4_units_type;

typedef bub::make_matrix<
    boost::fusion::vector<
        boost::fusion::vector<time_, length, velocity, dimensionless>,
        boost::fusion::vector<length, time_, velocity, dimensionless>,
        boost::fusion::vector<time_, velocity, length, dimensionless>,
        boost::fusion::vector<time_, length, dimensionless, velocity>
    >
>::type matrix_4x4_units_type;

struct derived_from_matrix_4x1_type :
    matrix_4x1_type
{};

struct derived_from_matrix_1x4_type :
    matrix_1x4_type
{};

struct derived_from_matrix_4x4_type :
    matrix_4x4_type
{};

namespace boost { namespace units_blas {
} }

int test_main (int, char *[])
{
    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::size<matrix_4x1_type>, boost::mpl::size_t<4> >));
    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::size<matrix_1x4_type>, boost::mpl::size_t<4> >));
    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::size<matrix_4x4_type>, boost::mpl::size_t<16> >));

    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::rows<matrix_4x1_units_type>, boost::mpl::size_t<4> >));
    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::rows<matrix_1x4_units_type>, boost::mpl::size_t<1> >));
    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::rows<matrix_4x4_units_type>, boost::mpl::size_t<4> >));

    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::columns<matrix_4x1_type>, boost::mpl::size_t<1> >));
    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::columns<matrix_1x4_type>, boost::mpl::size_t<4> >));
    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::columns<matrix_4x4_type>, boost::mpl::size_t<4> >));

    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::columns<matrix_4x1_units_type>, boost::mpl::size_t<1> >));
    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::columns<matrix_1x4_units_type>, boost::mpl::size_t<4> >));
    BOOST_MPL_ASSERT((boost::mpl::equal_to<bub::columns<matrix_4x4_units_type>, boost::mpl::size_t<4> >));

    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x1_type>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_1x4_type>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x4_type>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x1_type const>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x1_type volatile>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x1_type const volatile>));
    BOOST_MPL_ASSERT_NOT((bub::is_matrix<int>));
    BOOST_MPL_ASSERT_NOT((bub::is_matrix<double>));
    BOOST_MPL_ASSERT_NOT((bub::is_matrix<boost::fusion::vector<boost::fusion::vector<int, long, float, double> > >));

    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<matrix_4x1_type>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<matrix_1x4_type>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<matrix_4x4_type>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<derived_from_matrix_4x1_type>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<derived_from_matrix_1x4_type>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<derived_from_matrix_4x4_type>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<matrix_4x1_type const>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<matrix_4x1_type volatile>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<matrix_4x1_type const volatile>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<derived_from_matrix_4x1_type const>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<derived_from_matrix_4x1_type volatile>));
    BOOST_MPL_ASSERT((bub::is_or_is_derived_from_matrix<derived_from_matrix_4x1_type const volatile>));
    BOOST_MPL_ASSERT_NOT((bub::is_or_is_derived_from_matrix<int>));
    BOOST_MPL_ASSERT_NOT((bub::is_or_is_derived_from_matrix<double>));

    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x1_units_type>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_1x4_units_type>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x4_units_type>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x1_units_type const>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x1_units_type volatile>));
    BOOST_MPL_ASSERT((bub::is_matrix<matrix_4x1_units_type const volatile>));
    BOOST_MPL_ASSERT_NOT((bub::is_matrix<int>));
    BOOST_MPL_ASSERT_NOT((bub::is_matrix<double>));
    BOOST_MPL_ASSERT_NOT((bub::is_matrix<boost::fusion::vector<boost::fusion::vector<int, long, float, double> > >));

    BOOST_MPL_ASSERT_NOT((bub::is_square_matrix<matrix_4x1_type>));
    BOOST_MPL_ASSERT_NOT((bub::is_square_matrix<matrix_1x4_type>));
    BOOST_MPL_ASSERT((bub::is_square_matrix<matrix_4x4_type>));

    BOOST_MPL_ASSERT_NOT((bub::is_square_matrix<matrix_4x1_units_type>));
    BOOST_MPL_ASSERT_NOT((bub::is_square_matrix<matrix_1x4_units_type>));
    BOOST_MPL_ASSERT((bub::is_square_matrix<matrix_4x4_units_type>));

    BOOST_MPL_ASSERT((bub::is_vector<matrix_4x1_type>));
    BOOST_MPL_ASSERT_NOT((bub::is_vector<matrix_1x4_type>));
    BOOST_MPL_ASSERT_NOT((bub::is_vector<matrix_4x4_type>));

    BOOST_MPL_ASSERT((bub::is_vector<matrix_4x1_units_type>));
    BOOST_MPL_ASSERT_NOT((bub::is_vector<matrix_1x4_units_type>));
    BOOST_MPL_ASSERT_NOT((bub::is_vector<matrix_4x4_units_type>));

    BOOST_MPL_ASSERT_NOT((bub::is_transpose_vector<matrix_4x1_type>));
    BOOST_MPL_ASSERT((bub::is_transpose_vector<matrix_1x4_type>));
    BOOST_MPL_ASSERT_NOT((bub::is_transpose_vector<matrix_4x4_type>));

    BOOST_MPL_ASSERT_NOT((bub::is_transpose_vector<matrix_4x1_units_type>));
    BOOST_MPL_ASSERT((bub::is_transpose_vector<matrix_1x4_units_type>));
    BOOST_MPL_ASSERT_NOT((bub::is_transpose_vector<matrix_4x4_units_type>));

    return 0;
}
