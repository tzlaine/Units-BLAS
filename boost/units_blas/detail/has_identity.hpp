// boost.units_blas
//
// Copyright (C) 2010 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_HAS_IDENTITY_HPP
#define BOOST_UNITS_BLAS_HAS_IDENTITY_HPP

#include <boost/units_blas/traits.hpp>


namespace boost { namespace units_blas { namespace detail {

    // TODO: Implement in terms of row-wise DDV equality.
    template <typename Matrix>
    struct has_identity :
        is_matrix_or_derived<Matrix>::type
    {};

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_HAS_IDENTITY_HPP
