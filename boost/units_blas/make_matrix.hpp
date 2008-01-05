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
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/range_c.hpp>


namespace boost { namespace units_blas {

    template <typename Rows>
    struct canonical_matrix
    {
        typedef matrix<typename detail::deep_as_vector<Rows>::type> type;
    };

    namespace detail {

        template <typename T>
        struct uniform
        {
            template <typename Ignore>
            struct apply
            {
                typedef T type;
            };
        };

    } // namespace detail

    template <typename T, std::size_t Rows, std::size_t Columns>
    struct uniform_matrix
    {
        typedef typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                mpl::range_c<std::size_t, 0, Columns>,
                detail::uniform<T>
            >::type
        >::type row_type;

        typedef typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                mpl::range_c<std::size_t, 0, Rows>,
                detail::uniform<row_type>
            >::type
        >::type all_rows_type;

        typedef matrix<all_rows_type> type;
    };

    namespace detail {

        struct vectorize
        {
            template <class T>
            struct apply
            {
                typedef fusion::vector<T> type;
            };
        };

    } // namespace detail

    template <typename Elements>
    struct vector
    {
        typedef typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                Elements,
                detail::vectorize
            >::type
        >::type all_rows_type;

        typedef matrix<all_rows_type> type;
    };

    namespace detail {

        struct identity
        {
            template <class T>
            struct apply
            {
                typedef T type;
            };
        };

    } // namespace detail

    template <typename Elements>
    struct transpose_vector
    {
        typedef typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                Elements,
                detail::identity
            >::type
        >::type elements_as_fusion_vector;

        typedef matrix<fusion::vector<elements_as_fusion_vector> > type;
    };

    template <typename T, std::size_t N>
    struct uniform_vector
    {
        typedef typename uniform_matrix<T, N, 1>::type type;
    };

    template <typename T, std::size_t N>
    struct uniform_transpose_vector
    {
        typedef typename uniform_matrix<T, 1, N>::type type;
    };

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MAKE_MATRIX_HPP
