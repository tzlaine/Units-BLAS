// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_MAKE_MATRIX_HPP
#define BOOST_UNITS_BLAS_MAKE_MATRIX_HPP

#include <boost/units_blas/matrix.hpp>


namespace boost { namespace units_blas {

    namespace detail {

        template <typename T, std::size_t N>
        struct repeat_type
        {
            template <typename Seq>
            static auto call (Seq seq)
            { return repeat_type<T, N - 1>::call(push_back<T>(seq)); }
        };

        template <typename T>
        struct repeat_type<T, 0>
        {
            template <typename Seq>
            static auto call (Seq seq)
            { return seq; }
        };

    }

    /** Convenience metafunction that returns a @c matrix<> of dimension @c
        Rows x @c Columns, in which each element is of type @c T. */
    template <typename T, std::size_t Rows, std::size_t Columns>
    struct uniform_matrix_type
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        using type = decltype(
            detail::make_matrix<Rows, Columns>(
                detail::tuple_from_types(
                    detail::repeat_type<T, Rows * Columns>::call(
                        detail::type_sequence<>{}
                    )
                )
            )
        );
#else
        using type = detail::unspecified;
#endif
    };

    template <typename T, std::size_t Rows, std::size_t Columns>
    using uniform_matrix = typename uniform_matrix_type<T, Rows, Columns>::type;

    namespace detail {

        template <typename Tuple, std::size_t I, std::size_t N>
        struct vector_type_impl
        {
            template <typename Types>
            static constexpr auto call (Types)
            {
                using type = typename std::tuple_element<
                    I,
                    Tuple
                >::type;
                auto types = push_back<type>(Types{});
                return vector_type_impl<Tuple, I + 1, N>::call(types);
            }
        };

        template <typename Tuple, std::size_t N>
        struct vector_type_impl<Tuple, N, N>
        {
            template <typename ...T>
            struct impl
            {
                using type = matrix<T...>;
            };

            template <typename Seq>
            static constexpr auto call (Seq seq)
            { return typename impl<Seq>::type{}; }
        };

    }

    /** Convenience metafunction that returns a "vector" -- a @c matrix<> of
        dimension N x 1.  Specifically, the resulting @c matrix<>'s dimensions
        are size(@c Tuple) x 1, and each element (i, 0) is the ith type in
        @c Tuple. */
    template <typename Tuple>
    struct vector_type
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        using type =
            detail::vector_type_impl<Tuple, 0, std::tuple_size<Tuple>::value>;
#else
        using type = detail::unspecified;
#endif
    };

    template <typename Tuple>
    using vector = typename vector_type<Tuple>::type;

    /** Convenience metafunction that returns a "transpose vector" -- a @c
        matrix<> of dimension 1 x N.  Specifically, the resulting @c
        matrix<>'s dimensions are 1 x size(@c Tuple), and each element
        (0, i) is the ith type in @c Tuple. */
    template <typename Tuple>
    struct transpose_vector_type
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        using type =
            matrix_t<std::tuple<Tuple>, 1, std::tuple_size<Tuple>::value>;
#else
        using type = detail::unspecified;
#endif
    };

    template <typename Tuple>
    using transpose_vector = typename transpose_vector_type<Tuple>::type;

    /** Convenience metafunction that returns a "vector" -- a @c matrix<> of
        dimension N x 1, in which each element is of type @c T. */
    template <typename T, std::size_t N>
    struct uniform_vector_type
    {
        using type = uniform_matrix<T, N, 1>;
    };

    template <typename T, std::size_t N>
    using uniform_vector = typename uniform_vector_type<T, N>::type;

    /** Convenience metafunction that returns a "transpose vector" -- a @c
        matrix<> of dimension 1 x N, in which each element is of type @c T. */
    template <typename T, std::size_t N>
    struct uniform_transpose_vector_type
    {
        using type = uniform_matrix<T, 1, N>;
    };

    template <typename T, std::size_t N>
    using uniform_transpose_vector =
        typename uniform_transpose_vector_type<T, N>::type;

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MAKE_MATRIX_HPP
