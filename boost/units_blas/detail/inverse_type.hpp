// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_INVERSE_TYPE_HPP
#define BOOST_UNITS_BLAS_INVERSE_TYPE_HPP

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/detail/value_type.hpp>


namespace boost { namespace units_blas { namespace detail {

    // matrix inverse type
    template <typename Matrix, std::size_t I>
    struct inverse_transpose_element
    {
        using transpose_type =
            tuple_element_t<transpose_index<Matrix>(I), Matrix>;
        using type = decltype(
            std::declval<value_type_t<transpose_type>>() /
            std::declval<transpose_type>()
        );
    };

    template <typename Matrix, std::size_t I>
    using inverse_transpose_element_t =
        typename inverse_transpose_element<Matrix, I>::type;


    template <typename Matrix, std::size_t ...I>
    auto inverse_types (std::index_sequence<I...>)
    { return type_sequence<inverse_transpose_element_t<Matrix, I>...>{}; }

    template <typename Matrix>
    struct inverse_type
    {
        using type = decltype(
            make_matrix<Matrix::num_columns, Matrix::num_rows>(
                tuple_from_types(
                    inverse_types<Matrix>(
                        std::make_index_sequence<Matrix::num_elements>()
                    )
                )
            )
        );
    };

    template <typename Matrix>
    using inverse_type_t = typename inverse_type<Matrix>::type;

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_INVERSE_TYPE_HPP
