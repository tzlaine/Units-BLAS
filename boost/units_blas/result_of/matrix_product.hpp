// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_MATRIX_PRODUCT_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_MATRIX_PRODUCT_HPP

#include <boost/units_blas/result_of/transpose.hpp>
#include <boost/units_blas/result_of/detail/matrix_product.hpp>

#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/range_c.hpp>


namespace boost { namespace units_blas { namespace result_of {

    /** Returns the type that results from multiplying \a MatrixL and \a
        MatrixR.  \a MatrixL and \a MatrixR must be matrix<>s, and the number of
        columns in \a MatrixL must the same as the number of rows in \a MatrixR.
        Also, a matrix-product type must exist for \a MatrixL and \a MatrixR
        (some otherwise-suitable pairs of matrix<>s do not have a matrix-product
        that makes sense when their elements are unit types). */
    template <typename MatrixL, typename MatrixR>
    struct matrix_product
    {
        typedef matrix<
            typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, MatrixL::num_rows_t::value>,
                    detail::matrix_product_row<
                        typename MatrixL::value_types,
                        typename transpose<MatrixR>::type::value_types,
                        typename MatrixR::num_columns_t
                    >
                >::type
            >::type
        > type;
    };

    /** Returns the type that results from multiplying each element of \a Matrix
        by \a T.  \a Matrix must be a matrix<>. */
    template <typename Matrix, typename T>
    struct scalar_product
    {
        typedef matrix<
            typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, Matrix::num_rows_t::value>,
                    detail::scalar_product_row<Matrix, T, typename Matrix::num_columns_t>
                >::type
            >::type
        > type;
    };

    /** Returns the type that results from dividing each element of \a Matrix by
        \a T.  \a Matrix must be a matrix<>. */
    template <typename Matrix, typename T>
    struct scalar_quotient
    {
        typedef matrix<
            typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, Matrix::num_rows_t::value>,
                    detail::scalar_quotient_row<Matrix, T, typename Matrix::num_columns_t>
                >::type
            >::type
        > type;
    };

    /** Returns the type that results from adding \a MatrixL and \a MatrixR.  \a
        MatrixL and \a MatrixR must be matrix<>s with the same dimensions.
        Also, every sum MatrixL()(i, j) + MatrixR()(i, j) must be a valid
        operation. */
    template <typename MatrixL, typename MatrixR>
    struct matrix_element_sum
    {
        typedef matrix<
            typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, MatrixL::num_rows_t::value>,
                    detail::matrix_element_sum_row<
                        typename MatrixL::value_types,
                        typename MatrixR::value_types
                    >
                >::type
            >::type
        > type;
    };

    /** Returns the type that results from subtracting \a MatrixR from \a
        MatrixL.  \a MatrixL and \a MatrixR must be matrix<>s with the same
        dimensions.  Also, every difference MatrixL()(i, j) - MatrixR()(i, j)
        must be a valid operation. */
    template <typename MatrixL, typename MatrixR>
    struct matrix_element_difference
    {
        typedef matrix<
            typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, MatrixL::num_rows_t::value>,
                    detail::matrix_element_difference_row<
                        typename MatrixL::value_types,
                        typename MatrixR::value_types
                    >
                >::type
            >::type
        > type;
    };

    /** Returns the type that results from multiplying values in \a MatrixL by
        values in \a MatrixR, element-by-element.  \a MatrixL and \a MatrixR
        must be a matrix<>s with the same dimensions. */
    template <typename MatrixL, typename MatrixR>
    struct matrix_element_product
    {
        typedef matrix<
            typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, MatrixL::num_rows_t::value>,
                    detail::matrix_element_product_row<
                        typename MatrixL::value_types,
                        typename MatrixR::value_types
                    >
                >::type
            >::type
        > type;
    };

    /** Returns the type that results from dividing values in \a MatrixL by
        values in \a MatrixR, element-by-element.  \a MatrixL and \a MatrixR
        must be a matrix<>s with the same dimensions. */
    template <typename MatrixL, typename MatrixR>
    struct matrix_element_quotient
    {
        typedef matrix<
            typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, MatrixL::num_rows_t::value>,
                    detail::matrix_element_quotient_row<
                        typename MatrixL::value_types,
                        typename MatrixR::value_types
                    >
                >::type
            >::type
        > type;
    };

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_MATRIX_PRODUCT_HPP
