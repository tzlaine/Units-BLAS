// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_DETERMINANT_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_DETERMINANT_HPP

#include <boost/units_blas/result_of/slice.hpp>
#include <boost/units_blas/result_of/detail/matrix_product.hpp>

#include <boost/mpl/vector.hpp>


namespace boost { namespace units_blas { namespace result_of {

    template <typename Matrix, std::size_t Size = Matrix::num_rows_t::value>
    struct determinant;

    namespace detail {

        template <typename N, typename I>
        struct range_without_i
        {
            typedef typename mpl::range_c<std::size_t, 0, N::value> range;
            typedef typename mpl::fold<
                range,
                mpl::vector<>,
                mpl::if_<
                    mpl::equal_to<mpl::_2, I>,
                    mpl::_1,
                    mpl::push_back<mpl::_1, mpl::_2>
                >
            >::type type;
        };

        template <typename Matrix>
        struct top_row_element_subdeterminant
        {
            template <typename N>
            struct apply
            {
                typedef typename value_at<Matrix, size_t_<0>, N>::type top_row_element;
                typedef typename determinant<
                    typename slice<
                        Matrix,
                        typename range_without_i<
                            typename Matrix::num_rows_t,
                            size_t_<0>
                        >::type,
                        typename range_without_i<
                            typename Matrix::num_columns_t,
                            N
                        >::type
                    >::type
                >::type slice_determinant;
                typedef BOOST_TYPEOF((top_row_element() * slice_determinant())) type;
            };
        };

    } // namespace detail

#if BOOST_UNITS_BLAS_USE_INEXACT_DETERMINANT_TYPE

    template <typename Matrix>
    struct determinant
    {
        // This is not really correct.  This is the type that all terms in the
        // determinant expression must be operator-() / operator+()-compatible
        // with, but the result may be subtly different (it's value_type may
        // be float when it really should be double, for instance).  Consider
        // this case:
        // typedef boost::units_blas::matrix<
        //     boost::fusion::vector<
        //        boost::fusion::vector<float, double>,
        //        boost::fusion::vector<double, float>
        //     >
        // > matrix_t;
        // Note that determinant<matrix_t>::type will be float, because this
        // version of determinant<> uses only the diagonal elements.

        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));
        typedef typename mpl::transform_view<
            mpl::range_c<std::size_t, 0, Matrix::num_rows_t::value>,
            value_at<Matrix, mpl::_1, mpl::_1>
        >::type diagonal_types;

        typedef typename mpl::advance_c<
            typename mpl::begin<diagonal_types>::type,
            1
        >::type begin;
        typedef typename mpl::end<diagonal_types>::type end;
        typedef typename mpl::fold<
            mpl::iterator_range<begin, end>,
            typename mpl::front<diagonal_types>::type,
            detail::value_product
        >::type type;
    };

#else

    /** Returns the type that results from taking the determinant of @c
        Matrix.  @c Matrix must be a @c matrix<>, and must be square.  Also, a
        determinant type must exist for @c Matrix (some otherwise-suitable
        <c>matrix<>s</c> do not have a determinant that makes sense when their
        elements are unit types).  Note that the default implementation of
        this metafunction is (compile-time) O(n!) in the number of rows.  See
        the Configuration section of the documentation for a possible
        alternative. */
    template <typename Matrix, std::size_t Size>
    struct determinant
    {
#ifndef BOOST_UNITS_BLAS_DOXYGEN
        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));
        typedef typename mpl::transform_view<
            mpl::range_c<std::size_t, 0, Matrix::num_rows_t::value>,
            detail::top_row_element_subdeterminant<Matrix>
        >::type top_row_element_subdeterminants;

        typedef typename mpl::advance_c<
            typename mpl::begin<top_row_element_subdeterminants>::type,
            1
        >::type begin;
        typedef typename mpl::end<top_row_element_subdeterminants>::type end;
#endif

        typedef typename mpl::fold<
            mpl::iterator_range<begin, end>,
            typename mpl::front<top_row_element_subdeterminants>::type,
            typename mpl::lambda<detail::value_sum<mpl::_1, mpl::_2> >::type
        >::type type;
    };

#ifndef BOOST_UNITS_BLAS_DOXYGEN
    template <typename Matrix>
    struct determinant<Matrix, 1>
    {
        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));
        typedef typename value_at_c<Matrix, 0, 0>::type type;
    };

    template <typename Matrix>
    struct determinant<Matrix, 2>
    {
        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));
        typedef BOOST_TYPEOF((
            typename value_at_c<Matrix, 0, 0>::type() * typename value_at_c<Matrix, 1, 1>::type() -
            typename value_at_c<Matrix, 0, 1>::type() * typename value_at_c<Matrix, 1, 0>::type()
        )) type;
    };

    template <typename Matrix>
    struct determinant<Matrix, 3>
    {
        BOOST_MPL_ASSERT((mpl::equal_to<typename Matrix::num_rows_t,
                                        typename Matrix::num_columns_t>));
        typedef BOOST_TYPEOF((
            typename value_at_c<Matrix, 0, 0>::type() *
            (typename value_at_c<Matrix, 1, 1>::type() * typename value_at_c<Matrix, 2, 2>::type() -
             typename value_at_c<Matrix, 1, 2>::type() * typename value_at_c<Matrix, 2, 1>::type()) -
            typename value_at_c<Matrix, 0, 1>::type() *
            (typename value_at_c<Matrix, 1, 0>::type() * typename value_at_c<Matrix, 2, 2>::type() -
             typename value_at_c<Matrix, 1, 2>::type() * typename value_at_c<Matrix, 2, 0>::type()) +
            typename value_at_c<Matrix, 0, 2>::type() *
            (typename value_at_c<Matrix, 1, 0>::type() * typename value_at_c<Matrix, 2, 1>::type() -
             typename value_at_c<Matrix, 1, 1>::type() * typename value_at_c<Matrix, 2, 0>::type())
        )) type;
    };
#endif

#endif

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_DETERMINANT_HPP
