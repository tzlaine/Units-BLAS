// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_INVERSE_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_INVERSE_HPP

#include <boost/units_blas/result_of/transpose.hpp>
#include <boost/units_blas/result_of/detail/matrix_product.hpp>


namespace boost { namespace units_blas { namespace result_of {

    namespace detail {

        template <typename Matrix, typename Row, typename Column>
        struct inverse_element :
            mpl::lambda<
                value_inverse<
                    typename result_of::value_at<Matrix, Row, Column>::type
                >
            >::type
        {};

        template <typename Matrix, typename NumColumns>
        struct inverse_row
        {
            template <typename Row>
            struct apply
            {
                typedef typename fusion::result_of::as_vector<
                    typename mpl::transform_view<
                        mpl::range_c<std::size_t, 0, NumColumns::value>,
                        detail::inverse_element<Matrix, Row, mpl::_1>
                    >::type
                >::type type;
            };
        };

    } // namespace detail

    /** Returns the type of inverting @c Matrix.  @c Matrix must be a @c
        matrix<>, and must be square.  */
    template <typename Matrix>
    struct inverse
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));
#endif

        typedef matrix<
            typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, Matrix::num_rows_t::value>,
                    detail::inverse_row<
                        typename transpose<Matrix>::type,
                        typename Matrix::num_columns_t
                    >
                >::type
            >::type
        > type;
    };

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_INVERSE_HPP
