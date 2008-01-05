// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_SUM_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_SUM_HPP

#include <boost/units_blas/result_of/transpose.hpp>
#include <boost/units_blas/result_of/detail/matrix_product.hpp>

#include <boost/typeof/typeof.hpp>
#include <boost/mpl/assert.hpp>


namespace boost { namespace units_blas { namespace result_of {

    template <typename Matrix>
    struct sum
    {
        BOOST_MPL_ASSERT((
            mpl::or_<
                mpl::equal_to<typename Matrix::num_rows_t, size_t_<1> >,
                mpl::equal_to<typename Matrix::num_columns_t, size_t_<1> >
            >
        ));

        typedef typename mpl::eval_if<
            mpl::equal_to<typename Matrix::num_columns_t, size_t_<1> >,
            fusion::result_of::value_at_c<typename result_of::transpose<Matrix>::type::value_types, 0>,
            fusion::result_of::value_at_c<typename Matrix::value_types, 0>
        >::type element_types;

        typedef typename mpl::advance_c<
            typename mpl::begin<element_types>::type,
            1
        >::type begin;
        typedef typename mpl::end<element_types>::type end;
        typedef typename mpl::fold<
            mpl::iterator_range<begin, end>,
            typename mpl::front<element_types>::type,
            typename mpl::lambda<detail::value_sum<mpl::_1, mpl::_2> >::type
        >::type type;
    };

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_SUM_HPP
