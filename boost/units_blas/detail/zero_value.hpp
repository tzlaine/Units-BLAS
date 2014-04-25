// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DETAIL_ZERO_VALUE_HPP
#define BOOST_UNITS_BLAS_DETAIL_ZERO_VALUE_HPP

#include <boost/units/quantity.hpp>


namespace boost { namespace units_blas { namespace detail {

    template <typename T>
    struct zero_value
    {
        static T value ()
        { return T{0}; }
    };

    template <typename Unit, typename ValueType>
    struct zero_value<units::quantity<Unit, ValueType>>
    {
        static units::quantity<Unit, ValueType> value ()
        {
            return units::quantity<Unit, ValueType>::from_value(
                static_cast<ValueType>(0)
            );
        }
    };

} } }

#endif // BOOST_UNITS_BLAS_DETAIL_ZERO_VALUE_HPP
