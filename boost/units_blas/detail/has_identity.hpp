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

    template <typename T, typename Head, typename ...Tail>
    constexpr auto
    all_same_impl (std::pair<bool, type_sequence<Head, Tail...>> p)
    {
        return all_same_impl<T>(
            std::make_pair(
                p.first && std::is_same<T, Head>::value,
                type_sequence<Tail...>{}
            )
        );
    }

    template <typename T>
    constexpr auto all_same_impl (std::pair<bool, type_sequence<>> p)
    { return p.first; }

    template <typename Head, typename ...Tail>
    constexpr auto all_same (type_sequence<Head, Tail...> seq)
    { return all_same_impl<Head>(std::make_pair(true, seq)); }

    template <typename Matrix, std::size_t Row>
    struct row_ddv
    {
        static constexpr std::size_t origin = Row * Matrix::num_columns;

        template <std::size_t ...I>
        struct apply
        {
            using first = tuple_element_t<origin, Matrix>;
            using type = type_sequence<
                dimension_of_t<
                    decltype(
                        std::declval<first>() /
                        std::declval<tuple_element_t<origin + I, Matrix>>()
                    )
                >...
            >;
        };
    };

    template <typename Matrix, std::size_t ...I>
    constexpr auto row_ddvs (std::index_sequence<I...>)
    {
        return type_sequence<
            typename row_ddv<Matrix, I>::template apply<I...>::type...
        >{};
    }

    template <typename Matrix>
    struct matrix_has_identity
    {
        static constexpr auto all_ddvs =
            row_ddvs<Matrix>(std::make_index_sequence<Matrix::num_rows>());
        using type = std::integral_constant<bool, all_same(all_ddvs)>;
    };

    template <typename T>
    struct type_wrapper
    {
        using type = T;
    };

    template <typename Matrix>
    struct has_identity
    {
        using selector = std::conditional_t<
            is_matrix_or_derived<Matrix>::value,
            matrix_has_identity<Matrix>,
            type_wrapper<std::false_type>
        >;
        using type = typename selector::type;
        static const bool value = type::value;
    };

    template <typename Matrix>
    using has_identity_t = typename has_identity<Matrix>::type;

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_HAS_IDENTITY_HPP
