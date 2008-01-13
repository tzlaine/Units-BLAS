// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_TRANSPOSE_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_TRANSPOSE_HPP

#include <boost/units_blas/matrix_fwd.hpp>

//#include <boost/fusion/include/zip_view.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/range_c.hpp>


namespace boost { namespace units_blas { namespace result_of {

    namespace detail {

        template <typename ValueTypes, typename Row>
        struct transpose_element
        {
            template <typename Column>
            struct apply
            {
                typedef typename fusion::result_of::value_at<
                    typename fusion::result_of::value_at<ValueTypes, Column>::type,
                    Row
                >::type type;
            };
        };

        template <typename ValueTypes, typename NumColumns>
        struct transpose_row
        {
            template <typename Row>
            struct apply
            {
                typedef typename fusion::result_of::as_vector<
                    typename mpl::transform_view<
                        mpl::range_c<std::size_t, 0, NumColumns::value>,
                        transpose_element<ValueTypes, Row>
                    >::type
                >::type type;
            };
        };

    } // namespace detail

    /** Returns the type of taking the transpose of \a Matrix.  \a Matrix must
        be a matrix<>. */
    template <typename Matrix>
    struct transpose
    {
#if 0
        // TODO: The implementation outside this #if-block is more idiomatic
        // wrt the rest of the result_of code.  Test the compile time
        // performance of both befor removing either implementation.
        typedef typename Matrix::value_types value_types;
        typedef typename mpl::transform<
            value_types,
            add_reference<mpl::_1>
        >::type value_type_refs;
        typedef typename fusion::result_of::as_vector<
            fusion::zip_view<value_type_refs>
        >::type transpose_value_types;
        typedef matrix<transpose_value_types> type;
#endif
        typedef matrix<
            typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, Matrix::num_columns_t::value>,
                    detail::transpose_row<
                        typename Matrix::value_types,
                        typename Matrix::num_rows_t
                    >
                >::type
            >::type
        > type;
    };

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_TRANSPOSE_HPP
