// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_IDENTIY_MATRIX_HPP
#define BOOST_UNITS_BLAS_IDENTIY_MATRIX_HPP

#include <boost/units_blas/iterate.hpp>
#include <boost/units_blas/result_of/inverse.hpp>
#include <boost/units_blas/result_of/matrix_product.hpp>
#include <boost/units_blas/detail/simple_iteration.hpp>

namespace boost { namespace units_blas {

    /** Returns the type that can be used as the identity matrix for @c
        Matrix.  Note that not every matrix type has an identiy matrix
        type. */
    template <typename Matrix>
    struct identity_matrix
    {
        typedef typename result_of::matrix_product<
            Matrix,
            typename result_of::inverse<Matrix>::type
        >::type type;
    };

    /** Returns a matrix<> I whose diagonal elements are 1, and whose
        nondiagonal elements are 0, and whose element types are such that for
        any @c Matrix m, type of the expression m * I is assignable to a @c
        Matrix.  Further, for all i,j in [0, @c Matrix::num_rows_t::value), (m
        * I).at<i, j>() is within epsilon of m.at<i, j>().  Note that not
        every matrix type has an identiy matrix type. */
    template <typename Matrix>
    typename enable_if<
        mpl::equal_to<
            typename Matrix::num_rows_t,
            typename Matrix::num_columns_t
        >,
        typename identity_matrix<Matrix>::type
    >::type
    make_identity_matrix ()
    {
        typedef typename identity_matrix<Matrix>::type return_type;
        return_type retval;
        typedef fusion::vector<return_type &> ops;
        iterate<typename return_type::num_rows_t>(
            ops(retval), detail::identity_assign()
        );
        return retval;
    }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MATRIX_HPP
