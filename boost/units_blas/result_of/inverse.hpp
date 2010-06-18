// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_INVERSE_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_INVERSE_HPP

#include <boost/units_blas/result_of/determinant.hpp>
#include <boost/units_blas/result_of/detail/matrix_product.hpp>

#include <boost/mpl/vector.hpp>


namespace boost { namespace units_blas { namespace result_of {

    template <typename Matrix, std::size_t Size = Matrix::num_rows_t::value>
    struct inverse;

    namespace detail {

        template <typename Matrix, typename MatrixDeterminant, typename Row, typename Column>
        struct inverse_element_impl
        {
            typedef typename detail::value_inverse<
                typename determinant<Matrix>::type
            >::type inv_det;

            typedef typename determinant<
                typename slice<
                    Matrix,
                    typename range_without_i<
                        typename Matrix::num_rows_t,
                        Row
                    >::type,
                    typename range_without_i<
                        typename Matrix::num_columns_t,
                        Column
                    >::type
                >::type
            >::type subdeterminant;

            typedef BOOST_TYPEOF((subdeterminant() * inv_det())) type;
        };

        template <typename Matrix, typename MatrixDeterminant, typename Row, typename Column>
        struct inverse_element :
            mpl::lambda<
                inverse_element_impl<Matrix, MatrixDeterminant, Row, Column>
            >::type
        {};

        template <typename Matrix, typename MatrixDeterminant, typename NumColumns>
        struct inverse_row
        {
            template <typename Row>
            struct apply
            {
                typedef typename fusion::result_of::as_vector<
                    typename mpl::transform_view<
                        mpl::range_c<std::size_t, 0, NumColumns::value>,
                        detail::inverse_element<Matrix, MatrixDeterminant, Row, mpl::_1>
                    >::type
                >::type type;
            };
        };

    } // namespace detail

    // TODO: Doxygen comment
    template <typename Matrix, std::size_t Size>
    struct inverse
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));

        typedef typename determinant<Matrix>::type matrix_determinant;
#endif

        typedef matrix<
            typename fusion::result_of::as_vector<
            typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, Matrix::num_rows_t::value>,
                    detail::inverse_row<
                        Matrix,
                        matrix_determinant,
                        typename Matrix::num_columns_t
                    >
                >::type
            >::type
        > type;

    };

#ifndef BOOST_UNITS_BLAS_DOXYGEN
    template <typename Matrix>
    struct inverse<Matrix, 1>
    {
        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));
        typedef matrix<
            fusion::vector<
                fusion::vector<
                    typename detail::value_inverse<
                        typename value_at_c<Matrix, 0, 0>::type
                    >::type
                >
            >
        > type;
    };

    template <typename Matrix>
    struct inverse<Matrix, 2>
    {
        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));
        typedef typename detail::value_inverse<
            typename determinant<Matrix>::type
        >::type inv_det;
        typedef matrix<
            fusion::vector<
                fusion::vector<
                    BOOST_TYPEOF((typename value_at_c<Matrix, 1, 1>::type() * inv_det())),
                    BOOST_TYPEOF((typename value_at_c<Matrix, 1, 0>::type() * inv_det()))
                >,
                fusion::vector<
                    BOOST_TYPEOF((typename value_at_c<Matrix, 0, 1>::type() * inv_det())),
                    BOOST_TYPEOF((typename value_at_c<Matrix, 0, 0>::type() * inv_det()))
                >
            >
        > type;
    };

    template <typename Matrix>
    struct inverse<Matrix, 3>
    {
        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));
        typedef typename detail::value_inverse<
            typename determinant<Matrix>::type
        >::type inv_det;
        typedef matrix<
            fusion::vector<
                fusion::vector<
                    BOOST_TYPEOF((
                        (typename value_at_c<Matrix, 1, 2>::type() *
                         typename value_at_c<Matrix, 1, 2>::type() +
                         typename value_at_c<Matrix, 1, 1>::type() *
                         typename value_at_c<Matrix, 2, 2>::type()) * inv_det()
                    )),
                    BOOST_TYPEOF((
                        (typename value_at_c<Matrix, 0, 2>::type() *
                         typename value_at_c<Matrix, 2, 1>::type() +
                         typename value_at_c<Matrix, 0, 1>::type() *
                         typename value_at_c<Matrix, 2, 2>::type()) * inv_det()
                    )),
                    BOOST_TYPEOF((
                        (typename value_at_c<Matrix, 0, 2>::type() *
                         typename value_at_c<Matrix, 1, 1>::type() +
                         typename value_at_c<Matrix, 0, 1>::type() *
                         typename value_at_c<Matrix, 1, 2>::type()) * inv_det()
                    ))
                >,
                fusion::vector<
                    BOOST_TYPEOF((
                        (typename value_at_c<Matrix, 1, 2>::type() *
                         typename value_at_c<Matrix, 2, 0>::type() +
                         typename value_at_c<Matrix, 1, 0>::type() *
                         typename value_at_c<Matrix, 2, 2>::type()) * inv_det()
                    )),
                    BOOST_TYPEOF((
                        (typename value_at_c<Matrix, 0, 2>::type() *
                         typename value_at_c<Matrix, 2, 0>::type() +
                         typename value_at_c<Matrix, 0, 0>::type() *
                         typename value_at_c<Matrix, 2, 2>::type()) * inv_det()
                    )),
                    BOOST_TYPEOF((
                        (typename value_at_c<Matrix, 0, 2>::type() *
                         typename value_at_c<Matrix, 1, 0>::type() +
                         typename value_at_c<Matrix, 0, 0>::type() *
                         typename value_at_c<Matrix, 1, 2>::type()) * inv_det()
                    ))
                >,
                fusion::vector<
                    BOOST_TYPEOF((
                        (typename value_at_c<Matrix, 1, 1>::type() *
                         typename value_at_c<Matrix, 2, 0>::type() +
                         typename value_at_c<Matrix, 1, 0>::type() *
                         typename value_at_c<Matrix, 2, 1>::type()) * inv_det()
                    )),
                    BOOST_TYPEOF((
                        (typename value_at_c<Matrix, 0, 1>::type() *
                         typename value_at_c<Matrix, 2, 0>::type() +
                         typename value_at_c<Matrix, 0, 0>::type() *
                         typename value_at_c<Matrix, 2, 1>::type()) * inv_det()
                    )),
                    BOOST_TYPEOF((
                        (typename value_at_c<Matrix, 0, 1>::type() *
                         typename value_at_c<Matrix, 1, 0>::type() +
                         typename value_at_c<Matrix, 0, 0>::type() *
                         typename value_at_c<Matrix, 1, 1>::type()) * inv_det()
                    ))
                >
            >
        > type;
    };
#endif

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_INVERSE_HPP
