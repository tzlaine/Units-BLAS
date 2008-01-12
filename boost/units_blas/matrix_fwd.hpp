// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_MATRIX_FWD_HPP
#define BOOST_UNITS_BLAS_MATRIX_FWD_HPP

#include <boost/mpl/integral_c.hpp>


namespace boost { namespace units_blas {

    /** An mpl::integral_c<> whose value_type is std::size_t. */
    template <std::size_t N>
    struct size_t_ : mpl::integral_c<std::size_t, N> {};

    template <typename Rows>
    class matrix;

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MATRIX_FWD_HPP
