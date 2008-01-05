// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DETAIL_GET_VALUE_TYPE_HPP
#define BOOST_UNITS_BLAS_DETAIL_GET_VALUE_TYPE_HPP

#include <boost/units_blas/detail/has_value_type.hpp>

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>


namespace boost { namespace units_blas { namespace detail {

    template <typename T>
    struct value_type
    {
        typedef typename T::value_type type;
    };

    template <typename T>
    struct get_value_type
    {
        typedef typename mpl::eval_if<
            has_value_type<T>,
            value_type<T>,
            mpl::identity<T>
        >::type type;
    };

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_DETAIL_GET_VALUE_TYPE_HPP
