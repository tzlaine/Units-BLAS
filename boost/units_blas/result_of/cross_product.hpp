// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_CROSS_PRODUCT_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_CROSS_PRODUCT_HPP

#include <boost/units_blas/result_of/value_at.hpp>

#include <boost/typeof/typeof.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/or.hpp>


namespace boost { namespace units_blas { namespace result_of {

    namespace detail {

        template <typename VectorL, typename VectorR>
        struct vector_cross_product
        {
            typedef BOOST_TYPEOF((
                typename value_at_c<VectorL, 1, 0>::type() *
                typename value_at_c<VectorR, 2, 0>::type() -
                typename value_at_c<VectorL, 2, 0>::type() *
                typename value_at_c<VectorR, 1, 0>::type()
            )) element0;

            typedef BOOST_TYPEOF((
                typename value_at_c<VectorL, 2, 0>::type() *
                typename value_at_c<VectorR, 0, 0>::type() -
                typename value_at_c<VectorL, 0, 0>::type() *
                typename value_at_c<VectorR, 2, 0>::type()
            )) element1;

            typedef BOOST_TYPEOF((
                typename value_at_c<VectorL, 0, 0>::type() *
                typename value_at_c<VectorR, 1, 0>::type() -
                typename value_at_c<VectorL, 1, 0>::type() *
                typename value_at_c<VectorR, 0, 0>::type()
            )) element2;

            typedef matrix<
                fusion::vector<
                    fusion::vector<element0>,
                    fusion::vector<element1>,
                    fusion::vector<element2>
                >
            > type;
        };

        template <typename VectorL, typename VectorR>
        struct transpose_vector_cross_product
        {
            typedef BOOST_TYPEOF((
                typename value_at_c<VectorL, 0, 1>::type() *
                typename value_at_c<VectorR, 0, 2>::type() -
                typename value_at_c<VectorL, 0, 2>::type() *
                typename value_at_c<VectorR, 0, 1>::type()
            )) element0;

            typedef BOOST_TYPEOF((
                typename value_at_c<VectorL, 0, 2>::type() *
                typename value_at_c<VectorR, 0, 0>::type() -
                typename value_at_c<VectorL, 0, 0>::type() *
                typename value_at_c<VectorR, 0, 2>::type()
            )) element1;

            typedef BOOST_TYPEOF((
                typename value_at_c<VectorL, 0, 0>::type() *
                typename value_at_c<VectorR, 0, 1>::type() -
                typename value_at_c<VectorL, 0, 1>::type() *
                typename value_at_c<VectorR, 0, 0>::type()
            )) element2;

            typedef matrix<
                fusion::vector<
                    fusion::vector<element0, element1, element2>
                >
            > type;
        };

    } // namespace detail

    template <typename MatrixL, typename MatrixR>
    struct cross_product
    {
        BOOST_MPL_ASSERT((mpl::equal_to<typename MatrixL::num_rows_t,
                                        typename MatrixR::num_rows_t>));
        BOOST_MPL_ASSERT((mpl::equal_to<typename MatrixL::num_columns_t,
                                        typename MatrixR::num_columns_t>));
        BOOST_MPL_ASSERT((
            mpl::or_<
                mpl::and_<
                    mpl::equal_to<typename MatrixL::num_rows_t, size_t_<3> >,
                    mpl::equal_to<typename MatrixL::num_columns_t, size_t_<1> >
                >,
                mpl::and_<
                    mpl::equal_to<typename MatrixL::num_rows_t, size_t_<1> >,
                    mpl::equal_to<typename MatrixL::num_columns_t, size_t_<3> >
                >
            >
        ));

        typedef typename mpl::eval_if<
            mpl::equal_to<typename MatrixL::num_rows_t, size_t_<3> >,
            detail::vector_cross_product<MatrixL, MatrixR>,
            detail::transpose_vector_cross_product<MatrixL, MatrixR>
        >::type type;
    };

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_CROSS_PRODUCT_HPP
