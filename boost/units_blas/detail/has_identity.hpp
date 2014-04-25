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
#include <boost/units_blas/detail/lu.hpp>


namespace boost { namespace units_blas { namespace detail {

#if 0
    template <typename Matrix, std::size_t Row>
    struct row_ddv
    {
        constexpr std::size_t origin = Row * Matrix::num_columns;

        template <std::size_t ...I>
        constexpr auto call (std::index_sequence<I...>)
        {
            using first = tuple_element_t<origin, Matrix>;
            return type_sequence<
                dimension_of_t<
                    std::declval<first>() /
                    std::declval<tuple_element_t<origin + I, Matrix>>()
                >...
            >{};
        }
    };
#endif

    // TODO: Implement in terms of row-wise DDV equality.
    template <typename Matrix>
    struct has_identity :
        is_matrix_or_derived<Matrix>::type
    {};

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_HAS_IDENTITY_HPP
