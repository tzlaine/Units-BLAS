// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/traits.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

#include <boost/test/minimal.hpp>


namespace bub = boost::units_blas;

typedef bub::matrix<
    std::tuple<int>,
    std::tuple<long>,
    std::tuple<float>,
    std::tuple<double>
> matrix_4x1_type;

typedef bub::matrix<
    std::tuple<int, long, float, double>
> matrix_1x4_type;

#if 0 // TODO: This works around a bug in clang-3.4.
typedef bub::matrix<
    std::tuple<int, long, float, double>,
    std::tuple<long, int, float, double>,
    std::tuple<int, float, long, double>,
    std::tuple<int, long, double, float>
> matrix_4x4_type;
#else
typedef bub::matrix_t<
    std::tuple<int, long, float, double,
               long, int, float, double,
               int, float, long, double,
               int, long, double, float>,
    4,
    4
> matrix_4x4_type;
#endif

typedef boost::units::quantity<boost::units::si::time> time_;
typedef boost::units::quantity<boost::units::si::length> length;
typedef boost::units::quantity<boost::units::si::velocity> velocity;
typedef boost::units::quantity<boost::units::si::dimensionless> dimensionless;

typedef bub::matrix<
    std::tuple<time_>,
    std::tuple<length>,
    std::tuple<velocity>,
    std::tuple<dimensionless>
> matrix_4x1_units_type;

typedef bub::matrix<
    std::tuple<time_, length, velocity, dimensionless>
> matrix_1x4_units_type;

typedef bub::matrix<
    std::tuple<time_, length, velocity, dimensionless>,
    std::tuple<length, time_, velocity, dimensionless>,
    std::tuple<time_, velocity, length, dimensionless>,
    std::tuple<time_, length, dimensionless, velocity>
> matrix_4x4_units_type;

struct derived_from_matrix_4x1_type :
    matrix_4x1_type
{};

struct derived_from_matrix_1x4_type :
    matrix_1x4_type
{};

struct derived_from_matrix_4x4_type :
    matrix_4x4_type
{};

int test_main (int, char *[])
{
    static_assert(bub::is_matrix<matrix_4x1_type>::value, "");
    static_assert(bub::is_matrix<matrix_1x4_type>::value, "");
    static_assert(bub::is_matrix<matrix_4x4_type>::value, "");
    static_assert(bub::is_matrix<matrix_4x1_type const>::value, "");
    static_assert(bub::is_matrix<matrix_4x1_type volatile>::value, "");
    static_assert(bub::is_matrix<matrix_4x1_type const volatile>::value, "");
    static_assert(!bub::is_matrix<int>::value, "");
    static_assert(!bub::is_matrix<double>::value, "");

    static_assert(bub::is_matrix_or_derived<matrix_4x1_type>::value, "");
    static_assert(bub::is_matrix_or_derived<matrix_1x4_type>::value, "");
    static_assert(bub::is_matrix_or_derived<matrix_4x4_type>::value, "");
    static_assert(bub::is_matrix_or_derived<derived_from_matrix_4x1_type>::value, "");
    static_assert(bub::is_matrix_or_derived<derived_from_matrix_1x4_type>::value, "");
    static_assert(bub::is_matrix_or_derived<derived_from_matrix_4x4_type>::value, "");
    static_assert(bub::is_matrix_or_derived<matrix_4x1_type const>::value, "");
    static_assert(bub::is_matrix_or_derived<matrix_4x1_type volatile>::value, "");
    static_assert(bub::is_matrix_or_derived<matrix_4x1_type const volatile>::value, "");
    static_assert(bub::is_matrix_or_derived<derived_from_matrix_4x1_type const>::value, "");
    static_assert(bub::is_matrix_or_derived<derived_from_matrix_4x1_type volatile>::value, "");
    static_assert(bub::is_matrix_or_derived<derived_from_matrix_4x1_type const volatile>::value, "");
    static_assert(!bub::is_matrix_or_derived<int>::value, "");
    static_assert(!bub::is_matrix_or_derived<double>::value, "");

    static_assert(bub::is_matrix<matrix_4x1_units_type>::value, "");
    static_assert(bub::is_matrix<matrix_1x4_units_type>::value, "");
    static_assert(bub::is_matrix<matrix_4x4_units_type>::value, "");
    static_assert(bub::is_matrix<matrix_4x1_units_type const>::value, "");
    static_assert(bub::is_matrix<matrix_4x1_units_type volatile>::value, "");
    static_assert(bub::is_matrix<matrix_4x1_units_type const volatile>::value, "");
    static_assert(!bub::is_matrix<int>::value, "");
    static_assert(!bub::is_matrix<double>::value, "");

    return 0;
}
