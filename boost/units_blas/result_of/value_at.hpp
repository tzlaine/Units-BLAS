// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_VALUE_AT_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_VALUE_AT_HPP

#include <boost/units_blas/matrix_fwd.hpp>

#include <boost/fusion/sequence/intrinsic/value_at.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/less.hpp>


namespace boost { namespace units_blas { namespace result_of {

    /** Returns the type of the element at row \a I, column \a J of Matrix.
        Matrix must be a matrix<>. */
    template <typename Matrix, typename I, typename J>
    struct value_at
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        BOOST_MPL_ASSERT((mpl::less<I, typename Matrix::num_rows_t>));
        BOOST_MPL_ASSERT((mpl::less<J, typename Matrix::num_columns_t>));
#endif

        typedef typename fusion::result_of::value_at<
            typename fusion::result_of::value_at<
                typename Matrix::value_types,
                I
            >::type,
            J
        >::type type;
    };

    /** Returns the type of the element at row \a I, column \a J of Matrix.
        Matrix must be a matrix<>. */
    template <typename Matrix, std::size_t I, std::size_t J>
    struct value_at_c
    {
        typedef typename value_at<Matrix, size_t_<I>, size_t_<J> >::type type;
    };

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_VALUE_AT_HPP
