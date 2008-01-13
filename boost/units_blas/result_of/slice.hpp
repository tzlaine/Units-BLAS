// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_SLICE_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_SLICE_HPP

#include <boost/fusion/include/as_vector.hpp>
#include <boost/mpl/transform_view.hpp>


namespace boost { namespace units_blas { namespace result_of {

    namespace detail {

        template <typename Matrix, typename Columns>
        struct slice_row
        {
            template <typename Row>
            struct apply
            {
                typedef typename fusion::result_of::as_vector<
                    typename mpl::transform_view<
                        Columns,
                        value_at<Matrix, Row, mpl::_1>
                    >::type
                >::type type;
            };
        };

    } // namespace detail

    /** Returns a matrix<> type whose elements are composed of the rows and
        columns of \a Matrix specified in \a Rows and \a Columns.  \a Matrix
        must be a matrix<>.  \a Rows and \a Columns must be type sequences
        containing integral constants; all integral constants in \a Rows and \a
        Columns must fall within the numbers of rows and columns in Matrix,
        respectively. */
    template <typename Matrix, typename Rows, typename Columns>
    struct slice
    {
        typedef matrix<
            typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    Rows,
                    detail::slice_row<Matrix, Columns>
                >::type
            >::type
        > type;
    };

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_SLICE_HPP
