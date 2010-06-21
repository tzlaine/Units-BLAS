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
#include <boost/units_blas/detail/is_assignable.hpp>


namespace boost { namespace units_blas { namespace detail {

    template <typename Matrix>
    struct has_inverse :
        is_assignable<
            typename result_of::matrix_product<
                Matrix,
                typename identity_matrix<Matrix>::type
            >::type,
            Matrix
        >
    {};

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_HAS_INVERSE_HPP
