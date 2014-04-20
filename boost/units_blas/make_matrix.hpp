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
    struct make_uniform_matrix
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

#if 0
    /** Convenience metafunction that returns a "vector" -- a @c matrix<> of
        dimension N x 1.  Specifically, the resulting @c matrix<>'s dimensions
        are size(@c Tuple) x 1, and each element (i, 0) is the ith type in
        @c Tuple. */
    template <typename Tuple>
    struct make_vector
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        typedef typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                Tuple,
                mpl::bind1<detail::make_fusion_vector, mpl::_1>
            >::type
        >::type all_rows_type;

        typedef matrix<all_rows_type> type;
#else
        using type = detail::unspecified;
#endif
    };
#endif

    /** Convenience metafunction that returns a "transpose vector" -- a @c
        matrix<> of dimension 1 x N.  Specifically, the resulting @c
        matrix<>'s dimensions are 1 x size(@c Tuple), and each element
        (0, i) is the ith type in @c Tuple. */
    template <typename Tuple>
    struct make_transpose_vector
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        using type = matrix_t<std::tuple<Tuple>, 1, std::tuple_size<Tuple>::value>;
#else
        using type = detail::unspecified;
#endif
    };

    /** Convenience metafunction that returns a "vector" -- a @c matrix<> of
        dimension N x 1, in which each element is of type @c T. */
    template <typename T, std::size_t N>
    struct make_uniform_vector
    {
        using type = typename make_uniform_matrix<T, N, 1>::type;
    };

    /** Convenience metafunction that returns a "transpose vector" -- a @c
        matrix<> of dimension 1 x N, in which each element is of type @c T. */
    template <typename T, std::size_t N>
    struct make_uniform_transpose_vector
    {
        using type = typename make_uniform_matrix<T, 1, N>::type;
    };

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MAKE_MATRIX_HPP
