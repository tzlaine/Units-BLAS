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

    /** Convenience metafunction that returns a @c matrix<> of dimension @c
        Rows x @c Columns, in which each element is of type @c T. */
    template <typename T, std::size_t Rows, std::size_t Columns>
    struct uniform_matrix_type
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        using type = decltype(
            detail::make_matrix<Rows, Columns>(
                hana::repeat<hana::Tuple>(std::declval<T>(), hana::size_t<Rows * Columns>)
            )
        );
#else
        using type = detail::unspecified;
#endif
    };

    template <typename T, std::size_t Rows, std::size_t Columns>
    using uniform_matrix = typename uniform_matrix_type<T, Rows, Columns>::type;

    template <typename T, std::size_t Rows, std::size_t Columns>
    auto make_uniform_matrix ()
    { return uniform_matrix<T, Rows, Columns>{}; }

    /** Convenience metafunction that returns a "vector" -- a @c matrix<> of
        dimension N x 1.  Specifically, the resulting @c matrix<>'s dimensions
        are size(@c Tuple) x 1, and each element (i, 0) is the ith type in
        @c Tuple. */
    template <typename ...T>
    struct vector_type
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        using type = matrix_t<hana::_tuple<T...>, sizeof...(T), 1>;
#else
        using type = detail::unspecified;
#endif
    };

    template <typename ...T>
    using vector = typename vector_type<T...>::type;

    template <typename ...T>
    auto make_vector ()
    { return vector<T...>{}; }

    /** Convenience metafunction that returns a "transpose vector" -- a @c
        matrix<> of dimension 1 x N.  Specifically, the resulting @c
        matrix<>'s dimensions are 1 x size(@c Tuple), and each element
        (0, i) is the ith type in @c Tuple. */
    template <typename ...T>
    struct transpose_vector_type
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        using type = matrix<hana::_tuple<T...>>;
#else
        using type = detail::unspecified;
#endif
    };

    template <typename ...T>
    using transpose_vector = typename transpose_vector_type<T...>::type;

    template <typename ...T>
    auto make_transpose_vector ()
    { return transpose_vector<T...>{}; }

    /** Convenience metafunction that returns a "vector" -- a @c matrix<> of
        dimension N x 1, in which each element is of type @c T. */
    template <typename T, std::size_t N>
    struct uniform_vector_type
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        using type = uniform_matrix<T, N, 1>;
#else
        using type = detail::unspecified;
#endif
    };

    template <typename T, std::size_t N>
    using uniform_vector = typename uniform_vector_type<T, N>::type;

    template <typename T, std::size_t N>
    auto make_uniform_vector ()
    { return uniform_vector<T, N>{}; }

    /** Convenience metafunction that returns a "transpose vector" -- a @c
        matrix<> of dimension 1 x N, in which each element is of type @c T. */
    template <typename T, std::size_t N>
    struct uniform_transpose_vector_type
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        using type = uniform_matrix<T, 1, N>;
#else
        using type = detail::unspecified;
#endif
    };

    template <typename T, std::size_t N>
    using uniform_transpose_vector =
        typename uniform_transpose_vector_type<T, N>::type;

    template <typename T, std::size_t N>
    auto make_uniform_transpose_vector ()
    { return uniform_transpose_vector<T, N>{}; }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MAKE_MATRIX_HPP
