// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DETAIL_VALUE_TYPE_HPP
#define BOOST_UNITS_BLAS_DETAIL_VALUE_TYPE_HPP

#include <type_traits>
#include <utility>


namespace boost { namespace units_blas { namespace detail {

    auto has_value_type (...)
    { return std::false_type{}; }

    template <typename T>
    auto has_value_type (T, typename T::value_type* = 0)
    { return std::true_type{}; }

    template <typename T>
    struct value_type_wrapper
    {
        using value_type = T;
    };

    template <typename T>
    struct value_type
    {
        using selector = typename std::conditional<
            std::is_same<
                std::true_type,
                decltype(has_value_type(std::declval<T>()))
            >::value,
            T,
            value_type_wrapper<T>
        >::type;
        using type = typename selector::value_type;
    };

    template <typename T>
    using value_type_t = typename value_type<T>::type;

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_DETAIL_VALUE_TYPE_HPP
