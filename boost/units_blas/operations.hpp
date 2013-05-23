// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_OPERATIONS_HPP
#define BOOST_UNITS_BLAS_OPERATIONS_HPP

#include <boost/units_blas/config.hpp>
#include <boost/units_blas/result_of.hpp>
#include <boost/units_blas/traits.hpp>
#include <boost/units_blas/detail/get_value_type.hpp>
#include <boost/units_blas/detail/has_inverse.hpp>
#include <boost/units_blas/detail/iteration.hpp>
#include <boost/units_blas/detail/lu.hpp>
#include <boost/units_blas/detail/one_value.hpp>
#include <boost/units_blas/detail/simple_iteration.hpp>
#include <boost/units_blas/detail/zero_value.hpp>

#include <boost/array.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/units/cmath.hpp>

#include <cmath>


namespace boost { namespace units_blas {

    /** Returns a @c matrix<> consisting of only the rows and columns of @c m
        specified by @c Rows and @c Columns.  @c m must be a @c matrix<>.
        @c Rows and @c Columns must be type sequences containing integral
        constants; all integral constants in @c Rows and @c Columns must be
        less than the number of rows and columns in @c m, respectively.
        Note that duplication and order preservation are not checked for the
        constants in @c Rows and @c Columns.  It is therefore possible to use
        @c slice<>() to rearrange and/or duplicate rows and/or columns. */
    template <typename Rows, typename Columns, typename T>
    typename result_of::slice<matrix<T>, Rows, Columns>::type
    slice (matrix<T> const & m)
    {
        BOOST_MPL_ASSERT((mpl::less<mpl::int_<0>, mpl::size<Rows> >));
        BOOST_MPL_ASSERT((mpl::less<mpl::int_<0>, mpl::size<Columns> >));
        typedef typename result_of::slice<matrix<T>, Rows, Columns>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, matrix<T> const &> ops;
        iterate<size<result_type> >(
            ops(retval, m), detail::slice_assign<Rows, Columns>()
        );
        return retval;
    }

    /** Returns a const-preserved reference to the element at row @c I, column
        @c J of @c m.  @c m must be a @c matrix<>. */
    template <typename I, typename J, typename T>
    typename result_of::at<matrix<T>, I, J>::type
    at (matrix<T> & m)
    { return m.template at<I::value, J::value>(); }

    /** Returns a const-preserved reference to the element at row @c I, column
        @c J of @c m.  @c m must be a @c matrix<>. */
    template <std::size_t I, std::size_t J, typename T>
    typename result_of::at_c<matrix<T>, I, J>::type
    at_c (matrix<T> & m)
    { return m.template at<I, J>(); }

    /** Returns the tranpose of @c m.  @c m must be a @c matrix<>. */
    template <typename T>
    typename result_of::transpose<matrix<T> >::type
    transpose (matrix<T> const & m)
    {
        typedef typename result_of::transpose<matrix<T> >::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, matrix<T> const &> ops;
        iterate<size<matrix<T> > >(
            ops(retval, m), detail::transpose_assign()
        );
        return retval;
    }

    /** Returns the negation of @c m.  @c m must be a @c matrix<>. */
    template <typename T>
    matrix<T>
    neg (matrix<T> const & m)
    {
        matrix<T> retval;
        typedef fusion::vector<matrix<T> &, matrix<T> const &> ops;
        iterate<size<matrix<T> > >(
            ops(retval, m), detail::negate_assign()
        );
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the negation of @c m.  @c m must be a @c matrix<>. */
    template <typename T>
    matrix<T>
    operator- (matrix<T> const & m)
    { return neg(m); }

#endif

    /** Returns the elementwise sum of @c lhs and @c rhs. @c lhs and @c rhs
        must be <c>matrix<></c>s with the same dimensions.  Also, every sum
        <c>lhs(i, j) + rhs(i, j)</c> must be a valid operation.  */
    template <typename T, typename U>
    typename lazy_enable_if<
        is_same_shape_matrix<matrix<T>, matrix<U> >,
        result_of::matrix_element_sum<matrix<T>, matrix<U> >
    >::type
    sum (matrix<T> const & lhs, matrix<U> const & rhs)
    {
        typedef typename result_of::matrix_element_sum<matrix<T>, matrix<U> >::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, matrix<T> const &, matrix<U> const &> ops;
        iterate<size<matrix<T> > >(
            ops(retval, lhs, rhs), detail::matrix_matrix_elem_add_assign()
        );
        return retval;
    }

    /** Returns the elementwise difference of @c lhs and @c rhs. @c lhs and @c
        rhs must be <c>matrix<></c>s with the same dimensions.  Also, every
        difference <c>lhs(i, j) - rhs(i, j)</c> must be a valid operation.  */
    template <typename T, typename U>
    typename lazy_enable_if<
        is_same_shape_matrix<matrix<T>, matrix<U> >,
        result_of::matrix_element_sum<matrix<T>, matrix<U> >
    >::type
    diff (matrix<T> const & lhs, matrix<U> const & rhs)
    {
        typedef typename result_of::matrix_element_difference<matrix<T>, matrix<U> >::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, matrix<T> const &, matrix<U> const &> ops;
        iterate<size<matrix<T> > >(
            ops(retval, lhs, rhs), detail::matrix_matrix_elem_sub_assign()
        );
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the elementwise sum of @c lhs and @c rhs. @c lhs and @c rhs
        must be <c>matrix<></c>s with the same dimensions.  Also, every sum
        <c>lhs(i, j) + rhs(i, j)</c> must be a valid operation.  */
    template <typename T, typename U>
    typename lazy_enable_if<
        is_same_shape_matrix<matrix<T>, matrix<U> >,
        result_of::matrix_element_sum<matrix<T>, matrix<U> >
    >::type
    operator+ (matrix<T> const & lhs, matrix<U> const & rhs)
    { return sum(lhs, rhs); }

    /** Returns the elementwise difference of @c lhs and @c rhs. @c lhs and @c
        rhs must be <c>matrix<></c>s with the same dimensions.  Also, every
        difference <c>lhs(i, j) - rhs(i, j)</c> must be a valid operation.  */
    template <typename T, typename U>
    typename lazy_enable_if<
        is_same_shape_matrix<matrix<T>, matrix<U> >,
        result_of::matrix_element_difference<matrix<T>, matrix<U> >
    >::type
    operator- (matrix<T> const & lhs, matrix<U> const & rhs)
    { return diff(lhs, rhs); }

#endif

    /** Returns the matrix-product of @c lhs and @c rhs. @c lhs and @c rhs
        must be <c>matrix<></c>s, and the number of columns in @c lhs must the
        same as the number of rows in @c rhs.  Also, a matrix-product type
        must exist for @c lhs and @c rhs (some otherwise-suitable pairs of
        <c>matrix<></c>s do not have a matrix-product that makes sense when
        their elements are unit types). */
    template <typename T, typename U>
    typename lazy_enable_if<
        mpl::equal_to<
            typename matrix<T>::num_columns_t,
            typename matrix<U>::num_rows_t
        >,
        result_of::matrix_product<matrix<T>, matrix<U> >
    >::type
    prod (matrix<T> const & lhs, matrix<U> const & rhs)
    {
        typedef typename result_of::matrix_product<matrix<T>, matrix<U> >::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, matrix<T> const &, matrix<U> const &> ops;
        iterate<
            mpl::times<
                typename result_type::num_rows_t,
                typename result_type::num_columns_t,
                typename matrix<T>::num_columns_t
            >
        >(ops(retval, lhs, rhs), detail::matrix_matrix_mul_assign());
        return retval;
    }

    /** Returns the product of @c m and @c t.  @c m must be a @c matrix<>, and
        @c T must not be a @c matrix<>. */
    template <typename T, typename U>
    typename lazy_disable_if<
        is_or_is_derived_from_matrix<T>,
        result_of::scalar_product<matrix<U>, T>
    >::type
    prod (matrix<U> const & m, T const & t)
    {
        typedef typename result_of::scalar_product<matrix<U>, T>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, matrix<U> const &, T const &> ops;
        iterate<size<result_type> >(
            ops(retval, m, t), detail::matrix_scalar_mul_assign()
        );
        return retval;
    }

    /** Returns the product of @c t and @c m.  @c m must be a @c matrix<>, and
        @c T must not be a @c matrix<>. */
    template <typename T, typename U>
    typename lazy_disable_if<
        is_or_is_derived_from_matrix<T>,
        result_of::scalar_product<matrix<U>, T>
    >::type
    prod (T const & t, matrix<U> const & m)
    { return prod(m, t); }

    /** Returns the result of dividing @c m by @c t.  @c m must be a @c
        matrix<>, and @c T must not be a @c matrix<>. */
    template <typename T, typename U>
    typename lazy_disable_if<
        is_or_is_derived_from_matrix<T>,
        result_of::scalar_quotient<matrix<U>, T>
    >::type
    div (matrix<U> const & m, T const & t)
    {
        typedef typename result_of::scalar_quotient<matrix<U>, T>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, matrix<U> const &, T const &> ops;
        iterate<size<result_type> >(
            ops(retval, m, t), detail::matrix_scalar_div_assign()
        );
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the matrix-product of @c lhs and @c rhs. @c lhs and @c rhs
        must be <c>matrix<></c>s, and the number of columns in @c lhs must the
        same as the number of rows in @c rhs.  Also, a matrix-product type
        must exist for @c lhs and @c rhs (some otherwise-suitable pairs of
        <c>matrix<></c>s do not have a matrix-product that makes sense when
        their elements are unit types). */
    template <typename T, typename U>
    typename lazy_enable_if<
        mpl::equal_to<
            typename matrix<T>::num_columns_t,
            typename matrix<U>::num_rows_t
        >,
        result_of::matrix_product<matrix<T>, matrix<U> >
    >::type
    operator* (matrix<T> const & lhs, matrix<U> const & rhs)
    { return prod(lhs, rhs); }

    /** Returns the product of @c m and @c t.  @c m must be a @c matrix<>, and
        @c T must not be a @c matrix<>. */
    template <typename T, typename U>
    typename lazy_disable_if<
        is_or_is_derived_from_matrix<T>,
        result_of::scalar_product<matrix<U>, T>
    >::type
    operator* (matrix<U> const & m, T const & t)
    { return prod(m, t); }

    /** Returns the product of @c t and @c m.  @c m must be a @c matrix<>, and
        @c T must not be a @c matrix<>. */
    template <typename T, typename U>
    typename lazy_disable_if<
        is_or_is_derived_from_matrix<T>,
        result_of::scalar_product<matrix<U>, T>
    >::type
    operator* (T const & t, matrix<U> const & m)
    { return prod(m, t); }

    /** Returns the result of dividing @c m by @c t.  @c m must be a @c
        matrix<>, and @c T must not be a @c matrix<>. */
    template <typename T, typename U>
    typename lazy_disable_if<
        is_or_is_derived_from_matrix<T>,
        result_of::scalar_quotient<matrix<U>, T>
    >::type
    operator/ (matrix<U> const & m, T const & t)
    { return div(m, t); }

#endif

    /** Returns the elementwise multiplication of the elements of @c lhs by
        the elements of @c rhs.  @c lhs and @c rhs must be <c>matrix<></c>s
        with the same dimensions. */
    template <typename T, typename U>
    typename lazy_enable_if<
        is_same_shape_matrix<matrix<T>, matrix<U> >,
        result_of::matrix_element_product<matrix<T>, matrix<U> >
    >::type
    element_prod (matrix<T> const & lhs, matrix<U> const & rhs)
    {
        typedef typename result_of::matrix_element_product<matrix<T>, matrix<U> >::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, matrix<T> const &, matrix<U> const &> ops;
        iterate<size<result_type> >(
            ops(retval, lhs, rhs), detail::matrix_matrix_elem_mul_assign()
        );
        return retval;
    }

    /** Returns the elementwise division of the elements of @c lhs by the
        elements of @c rhs.  @c lhs and @c rhs must be <c>matrix<></c>s with
        the same dimensions. */
    template <typename T, typename U>
    typename lazy_enable_if<
        is_same_shape_matrix<matrix<T>, matrix<U> >,
        result_of::matrix_element_quotient<matrix<T>, matrix<U> >
    >::type
    element_div (matrix<T> const & lhs, matrix<U> const & rhs)
    {
        typedef typename result_of::matrix_element_quotient<matrix<T>, matrix<U> >::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, matrix<T> const &, matrix<U> const &> ops;
        iterate<size<result_type> >(
            ops(retval, lhs, rhs), detail::matrix_matrix_elem_div_assign()
        );
        return retval;
    }

    /** Swaps the values in @c lhs and @c rhs.  Note that this is an
        O(@c size<matrix<T> >::value) operation. */
    template <typename T>
    void
    swap (matrix<T> & lhs, matrix<T> & rhs)
    {
        typedef fusion::vector<matrix<T> &, matrix<T> &> ops;
        iterate<size<matrix<T> > >(
            ops(lhs, rhs), detail::swap()
        );
    }

    /** Returns the dot product of @c lhs and @c rhs.  @c VectorL and @c VectorR
        must be <c>matrix<></c>s, and must have the same dimensions.
        Additionally, both <c>matrix<></c>s must be "vectors" of length greater
        than 1. */
    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_vector<VectorL, VectorR>,
            mpl::not_equal_to<rows<VectorL>, mpl::size_t<1> >
        >,
        result_of::dot_product<VectorL, VectorR>
    >::type
    dot (VectorL const & lhs, VectorR const & rhs)
    {
        typedef typename result_of::detail::dot_product<
            typename fusion::result_of::value_at_c<
                typename result_of::transpose<VectorL>::type::value_types,
                0
            >::type,
            typename fusion::result_of::value_at_c<
                typename result_of::transpose<VectorR>::type::value_types,
                0
            >::type
        >::type result_type;
        result_type retval = detail::zero_value<result_type>::value();
        typedef fusion::vector<result_type &, VectorL const &, VectorR const &> ops;
        iterate<typename VectorL::num_rows_t>(
            ops(retval, lhs, rhs), detail::vector_vector_dot_product_assign()
        );
        return retval;
    }

    /** Returns the dot product of @c lhs and @c rhs.  @c VectorL and @c VectorR
        must be <c>matrix<></c>s, and must have the same dimensions.
        Additionally, both <c>matrix<></c>s must be "transpose vectors" of
        length greater than 1. */
    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_transpose_vector<VectorL, VectorR>,
            mpl::not_equal_to<columns<VectorL>, mpl::size_t<1> >
        >,
        result_of::dot_product<VectorL, VectorR>
    >::type
    dot (VectorL const & lhs, VectorR const & rhs)
    {
        typedef typename result_of::detail::dot_product<
            typename fusion::result_of::value_at_c<typename VectorL::value_types, 0>::type,
            typename fusion::result_of::value_at_c<typename VectorR::value_types, 0>::type
        >::type result_type;
        result_type retval = detail::zero_value<result_type>::value();
        typedef fusion::vector<result_type &, VectorL const &, VectorR const &> ops;
        iterate<typename VectorL::num_columns_t>(
            ops(retval, lhs, rhs), detail::transpose_vector_transpose_vector_dot_product_assign()
        );
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the dot product of @c lhs and @c rhs.  @c VectorL and @c VectorR
        must be <c>matrix<></c>s, and must have the same dimensions.
        Additionally, both <c>matrix<></c>s must be "vectors" of length greater
        than 1. */
    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_vector<VectorL, VectorR>,
            mpl::not_equal_to<rows<VectorL>, mpl::size_t<1> >
        >,
        result_of::dot_product<VectorL, VectorR>
    >::type
    operator* (VectorL const & lhs, VectorR const & rhs)
    { return dot(lhs, rhs); }

    /** Returns the dot product of @c lhs and @c rhs.  @c VectorL and @c VectorR
        must be <c>matrix<></c>s, and must have the same dimensions.
        Additionally, both <c>matrix<></c>s must be "transpose vectors" of
        length greater than 1. */
    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_transpose_vector<VectorL, VectorR>,
            mpl::not_equal_to<columns<VectorL>, mpl::size_t<1> >
        >,
        result_of::dot_product<VectorL, VectorR>
    >::type
    operator* (VectorL const & lhs, VectorR const & rhs)
    { return dot(lhs, rhs); }

#endif

    /** Returns the cross product of @c lhs with @c rhs.  @c VectorL and @c
        VectorR must both be <c>matrix<></c>s, and must both be 3 x 1
        "vectors".  Also, a cross product type must exist for @c VectorL and
        @c VectorR (some otherwise-suitable pairs of <c>matrix<></c>s do not
        have a cross product that makes sense when their elements are unit
        types). */
    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_vector<VectorL, VectorR>,
            mpl::equal_to<
                typename VectorL::num_rows_t,
                mpl::size_t<3>
            >
        >,
        result_of::cross_product<VectorL, VectorR>
    >::type
    cross (VectorL const & lhs, VectorR const & rhs)
    {
        typedef typename result_of::cross_product<VectorL, VectorR>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, VectorL const &, VectorR const &> ops;
        retval.template at<0, 0>() = lhs.template at<1, 0>() * rhs.template at<2, 0>() - lhs.template at<2, 0>() * rhs.template at<1, 0>();
        retval.template at<1, 0>() = lhs.template at<2, 0>() * rhs.template at<0, 0>() - lhs.template at<0, 0>() * rhs.template at<2, 0>();
        retval.template at<2, 0>() = lhs.template at<0, 0>() * rhs.template at<1, 0>() - lhs.template at<1, 0>() * rhs.template at<0, 0>();
        return retval;
    }

    /** Returns the cross product of @c lhs with @c rhs.  @c VectorL and @c
        VectorR must both be <c>matrix<></c>s, and must both be 1 x 3
        "transpose vectors".  Also, a cross product type must exist for @c
        VectorL and @c VectorR (some otherwise-suitable pairs of
        <c>matrix<></c>s do not have a cross product that makes sense when
        their elements are unit types). */
    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_transpose_vector<VectorL, VectorR>,
            mpl::equal_to<
                typename VectorL::num_columns_t,
                mpl::size_t<3>
            >
        >,
        result_of::cross_product<VectorL, VectorR>
    >::type
    cross (VectorL const & lhs, VectorR const & rhs)
    {
        typedef typename result_of::cross_product<VectorL, VectorR>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, VectorL const &, VectorR const &> ops;
        retval.template at<0, 0>() = lhs.template at<0, 1>() * rhs.template at<0, 2>() - lhs.template at<0, 2>() * rhs.template at<0, 1>();
        retval.template at<0, 1>() = lhs.template at<0, 2>() * rhs.template at<0, 0>() - lhs.template at<0, 0>() * rhs.template at<0, 2>();
        retval.template at<0, 2>() = lhs.template at<0, 0>() * rhs.template at<0, 1>() - lhs.template at<0, 1>() * rhs.template at<0, 0>();
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the cross product of @c lhs with @c rhs.  @c VectorL and @c
        VectorR must both be <c>matrix<></c>s, and must both be 3 x 1
        "vectors".  Also, a cross product type must exist for @c VectorL and
        @c VectorR (some otherwise-suitable pairs of <c>matrix<></c>s do not
        have a cross product that makes sense when their elements are unit
        types). */
    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_vector<VectorL, VectorR>,
            mpl::equal_to<
                typename VectorL::num_rows_t,
                mpl::size_t<3>
            >
        >,
        result_of::cross_product<VectorL, VectorR>
    >::type
    operator^ (VectorL const & lhs, VectorR const & rhs)
    { return cross(lhs, rhs); }

    /** Returns the cross product of @c lhs with @c rhs.  @c VectorL and @c
        VectorR must both be <c>matrix<></c>s, and must both be 1 x 3
        "transpose vectors".  Also, a cross product type must exist for @c
        VectorL and @c VectorR (some otherwise-suitable pairs of
        <c>matrix<></c>s do not have a cross product that makes sense when
        their elements are unit types). */
    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_transpose_vector<VectorL, VectorR>,
            mpl::equal_to<
                typename VectorL::num_columns_t,
                mpl::size_t<3>
            >
        >,
        result_of::cross_product<VectorL, VectorR>
    >::type
    operator^ (VectorL const & lhs, VectorR const & rhs)
    { return cross(lhs, rhs); }

#endif

    /** Returns the sum of all elements in @c v.  @c Vector must be a "vector"
        @c matrix<>.  Also, a sum type must exist for @c Vector (some
        otherwise-suitable <c>matrix<></c>s do not have a sum that makes sense
        when their elements are unit types). */
    template <typename Vector>
    typename lazy_enable_if<
        is_vector<Vector>,
        result_of::sum<Vector>
    >::type
    sum (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        result_type retval = detail::zero_value<result_type>::value();
        typedef fusion::vector<result_type &, Vector const &> ops;
        iterate<typename Vector::num_rows_t>(
            ops(retval, v), detail::vector_sum<false>()
        );
        return retval;
    }

    /** Returns the sum of all elements in @c v.  @c Vector must be a "transpose
        vector" @c matrix<> of length greater than 1.  Also, a sum type must
        exist for @c Vector (some otherwise-suitable <c>matrix<></c>s do not
        have a sum that makes sense when their elements are unit types). */
    template <typename Vector>
    typename lazy_enable_if<
        mpl::and_<
            is_transpose_vector<Vector>,
            mpl::not_equal_to<columns<Vector>, mpl::size_t<1> >
        >,
        result_of::sum<Vector>
    >::type
    sum (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        result_type retval = detail::zero_value<result_type>::value();
        typedef fusion::vector<result_type &, Vector const &> ops;
        iterate<typename Vector::num_columns_t>(
            ops(retval, v), detail::transpose_vector_sum<false>()
        );
        return retval;
    }

    /** Returns the sum of the absolute values of all elements in @c v.  @c
        Vector must be a "vector" @c matrix<>.  Also, a sum type must exist
        for @c Vector (some otherwise-suitable <c>matrix<></c>s do not have a
        sum that makes sense when their elements are unit types). */
    template <typename Vector>
    typename lazy_enable_if<
        is_vector<Vector>,
        result_of::sum<Vector>
    >::type
    norm_1 (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        result_type retval = detail::zero_value<result_type>::value();
        typedef fusion::vector<result_type &, Vector const &> ops;
        iterate<typename Vector::num_rows_t>(
            ops(retval, v), detail::vector_sum<true>()
        );
        return retval;
    }

    /** Returns the sum of the absolute values of all elements in @c v.  @c
        Vector must be a "transpose vector" @c matrix<> of length greater than
        1.  Also, a sum type must exist for @c Vector (some otherwise-suitable
        <c>matrix<></c>s do not have a sum that makes sense when their elements
        are unit types). */
    template <typename Vector>
    typename lazy_enable_if<
        mpl::and_<
            is_transpose_vector<Vector>,
            mpl::not_equal_to<columns<Vector>, mpl::size_t<1> >
        >,
        result_of::sum<Vector>
    >::type
    norm_1 (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        result_type retval = detail::zero_value<result_type>::value();
        typedef fusion::vector<result_type &, Vector const &> ops;
        iterate<typename Vector::num_columns_t>(
            ops(retval, v), detail::transpose_vector_sum<true>()
        );
        return retval;
    }

    /** Returns the square root of the sum of the squares of all elements in
        @c v.  @c Vector must be a "vector" @c matrix<>.  Also, a sum type
        must exist for @c Vector (some otherwise-suitable <c>matrix<></c>s do
        not have a sum that makes sense when their elements are unit
        types). */
    template <typename Vector>
    typename lazy_enable_if<
        is_vector<Vector>,
        result_of::sum<Vector>
    >::type
    norm_2 (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        typedef typename result_of::detail::value_product<result_type, result_type>::type temp_type;
        temp_type tmp = detail::zero_value<temp_type>::value();
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_rows_t>(
            ops(tmp, v), detail::vector_norm_2()
        );
        using std::sqrt;
        return sqrt(tmp);
    }

    /** Returns the square root of the sum of the squares of all elements in @c
        v.  @c Vector must be a "transpose vector" @c matrix<> of length greater
        than 1.  Also, a sum type must exist for @c Vector (some
        otherwise-suitable <c>matrix<></c>s do not have a sum that makes sense
        when their elements are unit types). */
    template <typename Vector>
    typename lazy_enable_if<
        mpl::and_<
            is_transpose_vector<Vector>,
            mpl::not_equal_to<columns<Vector>, mpl::size_t<1> >
        >,
        result_of::sum<Vector>
    >::type
    norm_2 (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        typedef typename result_of::detail::value_product<result_type, result_type>::type temp_type;
        temp_type tmp = detail::zero_value<temp_type>::value();
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_columns_t>(
            ops(tmp, v), detail::transpose_vector_norm_2()
        );
        using std::sqrt;
        return sqrt(tmp);
    }

    /** Returns the max of the absolute values of all elements in @c v.  @c
        Vector must be a "vector" @c matrix<>.  Also, @c operator< must must
        exist for all pairs of elements in @c Vector. */
    template <typename Vector>
    typename lazy_enable_if<
        is_vector<Vector>,
        result_of::sum<Vector>
    >::type
    norm_inf (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        typedef std::pair<result_type, std::size_t> temp_type;
        temp_type tmp(detail::abs_(static_cast<result_type>(v.template at<0, 0>())), 0);
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_rows_t>(
            ops(tmp, v), detail::vector_norm_inf()
        );
        return tmp.first;
    }

    /** Returns the max of the absolute values of all elements in @c v.  @c
        Vector must be a "transpose vector" @c matrix<> of length greater than
        1.  Also, @c operator< must must exist for all pairs of elements in @c
        Vector. */
    template <typename Vector>
    typename lazy_enable_if<
        mpl::and_<
            is_transpose_vector<Vector>,
            mpl::not_equal_to<columns<Vector>, mpl::size_t<1> >
        >,
        result_of::sum<Vector>
    >::type
    norm_inf (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        typedef std::pair<result_type, std::size_t> temp_type;
        temp_type tmp(detail::abs_(static_cast<result_type>(v.template at<0, 0>())), 0);
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_columns_t>(
            ops(tmp, v), detail::transpose_vector_norm_inf()
        );
        return tmp.first;
    }

    /** Returns the index of the first element in @c v equal to @c norm_inf(v).
        @c Vector must be a "vector" @c matrix<>.  Also, @c operator< must must
        exist for all pairs of elements in @c Vector. */
    template <typename Vector>
    typename enable_if<
        is_vector<Vector>,
        std::size_t
    >::type
    norm_inf_index (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type sum_type;
        typedef std::pair<sum_type, std::size_t> temp_type;
        temp_type tmp(detail::abs_(static_cast<sum_type>(v.template at<0, 0>())), 0);
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_rows_t>(
            ops(tmp, v), detail::vector_norm_inf()
        );
        return tmp.second;
    }

    /** Returns the index of the first element in @c v equal to @c norm_inf(v).
        @c Vector must be a "transpose vector" @c matrix<> of length greater
        than 1.  Also, @c operator< must must exist for all pairs of elements in
        @c Vector. */
    template <typename Vector>
    typename enable_if<
        mpl::and_<
            is_transpose_vector<Vector>,
            mpl::not_equal_to<columns<Vector>, mpl::size_t<1> >
        >,
        std::size_t
    >::type
    norm_inf_index (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type sum_type;
        typedef std::pair<sum_type, std::size_t> temp_type;
        temp_type tmp(detail::abs_(static_cast<sum_type>(v.template at<0, 0>())), 0);
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_columns_t>(
            ops(tmp, v), detail::transpose_vector_norm_inf()
        );
        return tmp.second;
    }

#ifndef BOOST_UNITS_BLAS_DOXYGEN
    template <typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_square_matrix<matrix<T> >,
            mpl::equal_to<
                typename matrix<T>::num_rows_t,
                mpl::size_t<1>
            >
        >,
        result_of::determinant<matrix<T> >
    >::type
    determinant (matrix<T> const & m)
    { return m.template at<0, 0>(); }

    template <typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_square_matrix<matrix<T> >,
            mpl::equal_to<
                typename matrix<T>::num_rows_t,
                mpl::size_t<2>
            >
        >,
        result_of::determinant<matrix<T> >
    >::type
    determinant (matrix<T> const & m)
    { return m.template at<0, 0>() * m.template at<1, 1>() - m.template at<0, 1>() * m.template at<1, 0>(); }

    template <typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_square_matrix<matrix<T> >,
            mpl::equal_to<
                typename matrix<T>::num_rows_t,
                mpl::size_t<3>
            >
        >,
        result_of::determinant<matrix<T> >
    >::type
    determinant (matrix<T> const & m)
    {
        return
            m.template at<0, 0>() * (m.template at<1, 1>() * m.template at<2, 2>() - m.template at<1, 2>() * m.template at<2, 1>()) -
            m.template at<0, 1>() * (m.template at<1, 0>() * m.template at<2, 2>() - m.template at<1, 2>() * m.template at<2, 0>()) +
            m.template at<0, 2>() * (m.template at<1, 0>() * m.template at<2, 1>() - m.template at<1, 1>() * m.template at<2, 0>());
    }
#endif

    /** Returns the determinant of @c m.  @c m must be a @c matrix<>, and must
        be square.  Also, a determinant type must exist for @c m (some
        otherwise-suitable <c>matrix<></c>s do not have a determinant that
        makes sense when their elements are unit types).  */
    template <typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_square_matrix<matrix<T> >,
            mpl::less<
                mpl::size_t<3>,
                typename matrix<T>::num_rows_t
            >
        >,
        result_of::determinant<matrix<T> >
    >::type
    determinant (matrix<T> const & m)
    {
        // If you're seeing an error here, you're trying to perform LU
        // decomposition on a matrix that has mixed units of the same
        // dimension (e.g. centimeters and meters in the same matrix).
        BOOST_MPL_ASSERT((detail::is_lu_decomposable<matrix<T> >));

        typedef typename result_of::determinant<matrix<T> >::type result_type;
        result_type retval = detail::one_value<result_type>::value();
        typedef typename detail::get_value_type<result_type>::type raw_value_type;
        typedef array<
            array<raw_value_type, matrix<T>::num_columns_t::value>,
            matrix<T>::num_rows_t::value
        > temp_matrix_type;
        temp_matrix_type temp_matrix;
        typedef fusion::vector<temp_matrix_type &, matrix<T> const &> ops1;
        iterate<size<matrix<T> > >(
            ops1(temp_matrix, m), detail::matrix_to_temp_assign<result_type>()
        );
        array<std::size_t, matrix<T>::num_rows_t::value> indices;
        try {
            retval *= detail::lu_decompose(temp_matrix, indices);
            typedef fusion::vector<result_type &, temp_matrix_type const &> ops2;
            iterate<typename matrix<T>::num_rows_t>(
                ops2(retval, temp_matrix), detail::accumulate_determinant()
            );
        } catch (singular_matrix const &) {
            retval = detail::zero_value<result_type>::value();
        }
        return retval;
    }

    /** Returns the inverse of @c m.  Throws @c singular_matrix if @c m is
        found to be singular.  @c m must be a @c matrix<>, and must be
        square. */
    template <typename T>
    typename lazy_enable_if<
        is_square_matrix<matrix<T> >,
        result_of::inverse<matrix<T> >
    >::type
    inverse (matrix<T> const & m)
    {
        // If you're seeing an error here, it's becaue you're trying to invert
        // a matrix for which no valid inverse type exists.  A valid inverse
        // of matrix type M must yield an identity matrix I such that M = M *
        // I is a valid operation.
        BOOST_MPL_ASSERT((detail::has_inverse<matrix<T> >));

        // If you're seeing an error here, you're trying to perform LU
        // decomposition on a matrix that has mixed units of the same
        // dimension (e.g. centimeters and meters in the same matrix).
        BOOST_MPL_ASSERT((detail::is_lu_decomposable<matrix<T> >));

        typedef typename result_of::inverse<matrix<T> >::type result_type;
        result_type retval;
        typedef typename result_of::determinant<matrix<T> >::type temp_value_type;
        typedef typename detail::get_value_type<temp_value_type>::type raw_value_type;
        typedef array<
            array<raw_value_type, matrix<T>::num_columns_t::value>,
            matrix<T>::num_rows_t::value
        > temp_matrix_type;
        temp_matrix_type temp_matrix;
        typedef fusion::vector<temp_matrix_type &, matrix<T> const &> ops1;
        iterate<size<matrix<T> > >(
            ops1(temp_matrix, m), detail::matrix_to_temp_assign<temp_value_type>()
        );
        typedef array<std::size_t, matrix<T>::num_rows_t::value> indices_type;
        indices_type indices;
        detail::lu_decompose(temp_matrix, indices);
        temp_matrix_type temp_result_matrix;
        typedef fusion::vector<temp_matrix_type &, temp_matrix_type const &, indices_type const &> ops2;
        iterate<typename matrix<T>::num_rows_t>(
            ops2(temp_result_matrix, temp_matrix, indices), detail::assign_inverted_column()
        );
        typedef fusion::vector<result_type &, temp_matrix_type const &> ops3;
        iterate<size<matrix<T> > >(
            ops3(retval, temp_result_matrix), detail::temp_to_matrix_assign()
        );
        return retval;
    }

    /** Returns the solution to the equation Ax = b in @c x.  Throws @c
        singular_matrix if @c A is found to be singular.  @c A must be a @c
        matrix<>, and must be square.  @c x and @c b must be "vector"
        <c>matrix<></c>s with the same dimensions, and must have a number of
        rows equal to the number of columns in @c A. */
    template <typename T, typename U, typename V>
    typename enable_if<
        mpl::and_<
            is_square_matrix<matrix<T> >,
            is_same_length_vector<matrix<U>, matrix<V> >,
            mpl::equal_to<
                typename matrix<T>::num_columns_t,
                typename matrix<U>::num_rows_t
            >
        >
    >::type
    solve (matrix<T> const & A, matrix<V> const & b, matrix<U> & x)
    {
        typedef BOOST_TYPEOF((matrix<T>() * matrix<U>())) a_times_b_type;
        BOOST_MPL_ASSERT((is_convertible<a_times_b_type, matrix<V> >));

        // If you're seeing an error here, you're trying to perform LU
        // decomposition on a matrix that has mixed units of the same
        // dimension (e.g. centimeters and meters in the same matrix).
        BOOST_MPL_ASSERT((detail::is_lu_decomposable<matrix<T> >));

        typedef typename result_of::determinant<matrix<T> >::type temp_value_type;
        typedef typename detail::get_value_type<temp_value_type>::type raw_value_type;
        typedef array<
            array<raw_value_type, matrix<T>::num_columns_t::value>,
            matrix<T>::num_rows_t::value
        > temp_matrix_type;
        temp_matrix_type temp_A;
        typedef fusion::vector<temp_matrix_type &, matrix<T> const &> ops1;
        iterate<size<matrix<T> > >(
            ops1(temp_A, A), detail::matrix_to_temp_assign<temp_value_type>()
        );
        array<std::size_t, matrix<T>::num_rows_t::value> indices;
        detail::lu_decompose(temp_A, indices);
        typedef array<raw_value_type, matrix<T>::num_rows_t::value> temp_vector_type;
        temp_vector_type temp_vector;
        typedef fusion::vector<temp_vector_type &, matrix<V> const &> ops2;
        iterate<typename matrix<T>::num_rows_t>(
            ops2(temp_vector, b), detail::matrix_to_temp_vector_assign<temp_value_type>()
        );
        detail::lu_substitute(temp_A, indices, temp_vector);
        typedef fusion::vector<matrix<U> &, temp_vector_type const &> ops3;
        iterate<typename matrix<U>::num_rows_t>(
            ops3(x, temp_vector), detail::temp_vector_to_matrix_assign()
        );
    }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_OPERATIONS_HPP
