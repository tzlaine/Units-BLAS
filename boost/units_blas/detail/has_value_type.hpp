// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DETAIL_HAS_VALUE_TYPE_HPP
#define BOOST_UNITS_BLAS_DETAIL_HAS_VALUE_TYPE_HPP

#include <boost/mpl/if.hpp>


namespace boost { namespace units_blas { namespace detail {

    typedef char has_value_type_yes;
    typedef char (&has_value_type_no)[2];

    template <class T>
    inline has_value_type_yes has_value_type_ (const T &, typename T::value_type * = 0);
    inline has_value_type_no has_value_type_ (...);

    template <typename T>
    struct has_value_type :
        mpl::if_c<
            sizeof(has_value_type_(T())) == sizeof(has_value_type_yes),
            mpl::true_,
            mpl::false_
        >::type
    {};

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_DETAIL_HAS_VALUE_TYPE_HPP
