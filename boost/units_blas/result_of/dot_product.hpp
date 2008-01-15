// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_DOT_PRODUCT_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_DOT_PRODUCT_HPP

#include <boost/units_blas/result_of/transpose.hpp>
#include <boost/units_blas/result_of/detail/matrix_product.hpp>


namespace boost { namespace units_blas { namespace result_of {

    /** Returns the type of the dot product of @c MatrixL and @c MatrixR.  @c
        MatrixL and @c MatrixR must be <c>matrix<></c>s, and must have the
        same dimensions.  Additionally, both <c>matrix<></c>s must be
        "vectors", or "transpose vectors". */
    template <typename MatrixL, typename MatrixR>
    struct dot_product
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        BOOST_MPL_ASSERT((mpl::equal_to<typename MatrixL::num_rows_t,
                                        typename MatrixR::num_rows_t>));
        BOOST_MPL_ASSERT((mpl::equal_to<typename MatrixL::num_columns_t,
                                        typename MatrixR::num_columns_t>));
        BOOST_MPL_ASSERT((
            mpl::or_<
                mpl::equal_to<typename MatrixL::num_rows_t, size_t_<1> >,
                mpl::equal_to<typename MatrixL::num_columns_t, size_t_<1> >
            >
        ));
#endif

        typedef typename mpl::eval_if<
            mpl::equal_to<typename MatrixL::num_columns_t, size_t_<1> >,
            result_of::detail::dot_product<
                typename fusion::result_of::value_at_c<
                    typename result_of::transpose<MatrixL>::type::value_types,
                    0
                >::type,
                typename fusion::result_of::value_at_c<
                    typename result_of::transpose<MatrixR>::type::value_types,
                    0
                >::type
            >,
            result_of::detail::dot_product<
                typename fusion::result_of::value_at_c<
                    typename MatrixL::value_types,
                    0
                >::type,
                typename fusion::result_of::value_at_c<
                    typename MatrixR::value_types,
                    0
                >::type
            >
        >::type type;
    };

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_DOT_PRODUCT_HPP
