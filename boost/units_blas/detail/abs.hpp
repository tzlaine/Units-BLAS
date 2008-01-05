// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DETAIL_ABS_HPP
#define BOOST_UNITS_BLAS_DETAIL_ABS_HPP

#include <boost/units_blas/detail/zero_value.hpp>


namespace boost { namespace units_blas { namespace detail {

    template <typename T>
    T abs_ (T const & t)
    { return (zero_value<T>::value() < t) ? t : -t; }

} } }  // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_DETAIL_ABS_HPP
