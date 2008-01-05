// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_MATRIX_FWD_HPP
#define BOOST_UNITS_BLAS_MATRIX_FWD_HPP

#include <boost/fusion/container/vector/limits.hpp>
#include <boost/preprocessor/repetition/enum_params_with_a_default.hpp>
#include <boost/mpl/integral_c.hpp>


namespace boost { namespace fusion {
    struct void_;
} }

namespace boost { namespace units_blas {

    template <std::size_t N>
    struct size_t_ : mpl::integral_c<std::size_t, N> {};

    template <
        BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(
            FUSION_MAX_VECTOR_SIZE,
            typename T,
            fusion::void_
        )
    >
    struct row;

    template <
        BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(
            FUSION_MAX_VECTOR_SIZE,
            typename T,
            fusion::void_
        )
    >
    struct all_rows;

    template <typename Rows>
    class matrix;

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MATRIX_FWD_HPP
