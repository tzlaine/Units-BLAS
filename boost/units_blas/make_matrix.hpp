// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_MAKE_MATRIX_HPP
#define BOOST_UNITS_BLAS_MAKE_MATRIX_HPP

#include <boost/units_blas/matrix_fwd.hpp>

#include <boost/fusion/include/as_vector.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/range_c.hpp>


namespace boost { namespace units_blas {

    /** Convenience metafunction that returns @c matrix<X>, where @c X is the
        template parameter @c Rows, converted to a @c fusion::vector of
        @c fusion::vectors. */
    template <typename Rows>
    struct make_matrix
    {
        typedef matrix<typename detail::deep_as_vector<Rows>::type> type;
    };

    /** Convenience metafunction that returns a @c matrix<> of dimension @c
        Rows x @c Columns, in which each element is of type @c T. */
    template <typename T, std::size_t Rows, std::size_t Columns>
    struct uniform_matrix
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        typedef typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                mpl::range_c<std::size_t, 0, Columns>,
                mpl::bind1<mpl::quote1<mpl::identity>, T>
            >::type
        >::type row_type;

        typedef typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                mpl::range_c<std::size_t, 0, Rows>,
                mpl::bind1<mpl::quote1<mpl::identity>, row_type>
            >::type
        >::type all_rows_type;

        typedef matrix<all_rows_type> type;
#else
        typedef detail::unspecified type;
#endif
    };

    namespace detail {
        struct make_fusion_vector
        {
            template <
                BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(
                    FUSION_MAX_VECTOR_SIZE, typename T, fusion::void_)
            >
            struct apply
            {
                typedef fusion::vector<BOOST_PP_ENUM_PARAMS(FUSION_MAX_VECTOR_SIZE, T)> type;
            };
        };
    }

    /** Convenience metafunction that returns a "vector" -- a @c matrix<> of
        dimension N x 1.  Specifically, the resulting @c matrix<>'s dimensions
        are size(@c Elements) x 1, and each element (i, 0) is the ith type in
        @c Elements. */
    template <typename Elements>
    struct vector
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        typedef typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                Elements,
                mpl::bind1<detail::make_fusion_vector, mpl::_1>
            >::type
        >::type all_rows_type;

        typedef matrix<all_rows_type> type;
#else
        typedef detail::unspecified type;
#endif
    };

    /** Convenience metafunction that returns a "transpose vector" -- a @c
        matrix<> of dimension 1 x N.  Specifically, the resulting @c
        matrix<>'s dimensions are 1 x size(@c Elements), and each element
        (0, i) is the ith type in @c Elements. */
    template <typename Elements>
    struct transpose_vector
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        typedef typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                Elements,
                mpl::quote1<mpl::identity>
            >::type
        >::type elements_as_fusion_vector;

        typedef matrix<fusion::vector<elements_as_fusion_vector> > type;
#else
        typedef detail::unspecified type;
#endif
    };

    /** Convenience metafunction that returns a "vector" -- a @c matrix<> of
        dimension N x 1, in which each element is of type @c T. */
    template <typename T, std::size_t N>
    struct uniform_vector
    {
        typedef typename uniform_matrix<T, N, 1>::type type;
    };

    /** Convenience metafunction that returns a "transpose vector" -- a @c
        matrix<> of dimension 1 x N, in which each element is of type @c T. */
    template <typename T, std::size_t N>
    struct uniform_transpose_vector
    {
        typedef typename uniform_matrix<T, 1, N>::type type;
    };

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MAKE_MATRIX_HPP
