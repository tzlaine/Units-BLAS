// boost.units_blas
//
// Copyright (C) 2010 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_HAS_INVERSE_HPP
#define BOOST_UNITS_BLAS_HAS_INVERSE_HPP

#include <boost/units_blas/identity_matrix.hpp>


namespace boost { namespace units_blas { namespace detail {

    template <typename Matrix, typename Enable = void>
    struct has_inverse :
        std::false_type
    {};

    template <typename Matrix>
    struct has_inverse<
        Matrix,
        decltype(
            std::declval<Matrix>() *
            std::declval<typename identity_matrix<Matrix>::type>() // TODO: Should be right identity.
        )
    > :
        std::true_type
    {};

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_HAS_INVERSE_HPP
