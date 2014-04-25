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
    template <typename Matrix, std::size_t I, std::size_t N>
    struct inverse_transpose_types_impl
    {
        template <typename Types>
        static constexpr auto call (Types types)
        {
            constexpr std::size_t row = I / Matrix::num_columns;
            constexpr std::size_t column = I % Matrix::num_columns;
            constexpr std::size_t transpose_i =
                column * Matrix::num_rows + row;
            using type = typename std::tuple_element<
                transpose_i,
                typename Matrix::value_types
            >::type;
            using inverse_type = decltype(
                std::declval<typename value_type<type>::type>() /
                std::declval<type>()
            );
            return inverse_transpose_types_impl<Matrix, I + 1, N>::call(
                push_back<inverse_type>(types)
            );
        }
    };

    template <typename Matrix, std::size_t N>
    struct inverse_transpose_types_impl<Matrix, N, N>
    {
        template <typename Seq>
        static constexpr auto call (Seq seq)
        { return seq; }
    };

    template <typename Matrix>
    struct inverse_type
    {
        using type = decltype(
            make_matrix<Matrix::num_columns, Matrix::num_rows>(
                tuple_from_types(
                    inverse_transpose_types_impl<
                        Matrix,
                        0,
                        Matrix::num_elements
                    >::call(type_sequence<>{})
                )
            )
        );
    };

    template <typename Matrix>
    using inverse_type_t = typename inverse_type<Matrix>::type;

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_INVERSE_TYPE_HPP
