// boost.units_blas
//
// Copyright (C) 2010 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_IS_ASSIGNABLE_HPP
#define BOOST_UNITS_BLAS_IS_ASSIGNABLE_HPP

#include <boost/units_blas/result_of/value_at.hpp>


namespace boost { namespace units_blas { namespace detail {

    template <typename MatrixL, typename MatrixR>
    struct is_assignable_impl
    {
        template <typename State, typename Index>
        struct apply
        {
            typedef typename result_of::value_at_c<
                MatrixL,
                Index::value / MatrixL::num_columns_t::value,
                Index::value % MatrixL::num_columns_t::value
            >::type left_type;
            typedef typename result_of::value_at_c<
                MatrixR,
                Index::value / MatrixL::num_columns_t::value,
                Index::value % MatrixL::num_columns_t::value
            >::type right_type;
            typedef typename mpl::and_<
                State,
                mpl::or_<
                    is_same<left_type, right_type>,
                    is_convertible<left_type, right_type>
                >
            >::type type;
        };
    };

    template <typename MatrixL, typename MatrixR>
    struct is_assignable :
        mpl::and_<
            mpl::equal_to<
                typename MatrixL::num_rows_t,
                typename MatrixR::num_rows_t
            >,
            mpl::equal_to<
                typename MatrixL::num_columns_t,
                typename MatrixR::num_columns_t
            >,
            mpl::fold<
                mpl::range_c<std::size_t, 0, MatrixL::num_rows_t::value * MatrixL::num_columns_t::value>,
                mpl::true_,
                is_assignable_impl<MatrixL, MatrixR>
            >
        >
    {};

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_IS_ASSIGNABLE_HPP
