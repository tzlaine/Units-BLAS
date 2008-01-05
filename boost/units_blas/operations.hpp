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
#include <boost/units_blas/detail/simple_iteration.hpp>
#include <boost/units_blas/detail/iteration.hpp>
#include <boost/units_blas/detail/lu.hpp>
#include <boost/units_blas/detail/one_value.hpp>
#include <boost/units_blas/detail/zero_value.hpp>

#include <boost/array.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/units/cmath.hpp>

#include <cmath>


namespace boost { namespace units_blas {

    template <typename Rows, typename Columns, typename Matrix>
    typename lazy_enable_if<
        is_matrix<Matrix>,
        result_of::slice<Matrix, Rows, Columns>
    >::type
    slice (Matrix const & m)
    {
        BOOST_MPL_ASSERT((mpl::less<mpl::int_<0>, mpl::size<Rows> >));
        BOOST_MPL_ASSERT((mpl::less<mpl::int_<0>, mpl::size<Columns> >));
        typedef typename result_of::slice<Matrix, Rows, Columns>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, Matrix const &> ops;
        iterate<size<result_type> >(
            ops(retval, m), detail::slice_assign<Rows, Columns>()
        );
        return retval;
    }

    template <typename I, typename J, typename Matrix>
    typename lazy_enable_if<
        is_matrix<Matrix>,
        result_of::at<Matrix, I, J>
    >::type
    at (Matrix & m)
    { return m.at<I::value, J::value>(); }

    template <std::size_t I, std::size_t J, typename Matrix>
    typename lazy_enable_if<
        is_matrix<Matrix>,
        result_of::at_c<Matrix, I, J>
    >::type
    at_c (Matrix & m)
    { return m.at<I, J>(); }

    template <typename Matrix>
    typename lazy_enable_if<
        is_matrix<Matrix>,
        result_of::transpose<Matrix>
    >::type
    transpose (Matrix const & m)
    {
        typedef typename result_of::transpose<Matrix>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, Matrix const &> ops;
        iterate<size<Matrix> >(
            ops(retval, m), detail::transpose_assign()
        );
        return retval;
    }

    template <typename Matrix>
    typename enable_if<
        is_matrix<Matrix>,
        Matrix
    >::type
    neg (Matrix const & m)
    {
        Matrix retval;
        typedef fusion::vector<Matrix &, Matrix const &> ops;
        iterate<size<Matrix> >(
            ops(retval, m), detail::negate_assign()
        );
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    template <typename Matrix>
    typename enable_if<
        is_matrix<Matrix>,
        Matrix
    >::type
    operator- (Matrix const & m)
    { return neg(m); }

#endif

    template <typename MatrixL, typename MatrixR>
    typename lazy_enable_if<
        is_same_shape_matrix<MatrixL, MatrixR>,
        result_of::matrix_element_sum<MatrixL, MatrixR>
    >::type
    sum (MatrixL const & lhs, MatrixR const & rhs)
    {
        typedef typename result_of::matrix_element_sum<MatrixL, MatrixR>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, MatrixL const &, MatrixR const &> ops;
        iterate<size<MatrixL> >(
            ops(retval, lhs, rhs), detail::matrix_matrix_elem_add_assign()
        );
        return retval;
    }

    template <typename MatrixL, typename MatrixR>
    typename lazy_enable_if<
        is_same_shape_matrix<MatrixL, MatrixR>,
        result_of::matrix_element_sum<MatrixL, MatrixR>
    >::type
    diff (MatrixL const & lhs, MatrixR const & rhs)
    {
        typedef typename result_of::matrix_element_difference<MatrixL, MatrixR>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, MatrixL const &, MatrixR const &> ops;
        iterate<size<MatrixL> >(
            ops(retval, lhs, rhs), detail::matrix_matrix_elem_sub_assign()
        );
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    template <typename MatrixL, typename MatrixR>
    typename lazy_enable_if<
        is_same_shape_matrix<MatrixL, MatrixR>,
        result_of::matrix_element_sum<MatrixL, MatrixR>
    >::type
    operator+ (MatrixL const & lhs, MatrixR const & rhs)
    { return sum(lhs, rhs); }

    template <typename MatrixL, typename MatrixR>
    typename lazy_enable_if<
        is_same_shape_matrix<MatrixL, MatrixR>,
        result_of::matrix_element_difference<MatrixL, MatrixR>
    >::type
    operator- (MatrixL const & lhs, MatrixR const & rhs)
    { return diff(lhs, rhs); }

#endif

    template <typename MatrixL, typename MatrixR>
    typename lazy_enable_if<
        mpl::and_<
            is_matrix<MatrixL>,
            is_matrix<MatrixR>,
            mpl::equal_to<
                typename MatrixL::num_columns_t,
                typename MatrixR::num_rows_t
            >
        >,
        result_of::matrix_product<MatrixL, MatrixR>
    >::type
    prod (MatrixL const & lhs, MatrixR const & rhs)
    {
        typedef typename result_of::matrix_product<MatrixL, MatrixR>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, MatrixL const &, MatrixR const &> ops;
        iterate<
            mpl::times<
                typename result_type::num_rows_t,
                typename result_type::num_columns_t,
                typename MatrixL::num_columns_t
            >
        >(ops(retval, lhs, rhs), detail::matrix_matrix_mul_assign());
        return retval;
    }

    template <typename Matrix, typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_matrix<Matrix>,
            mpl::not_<is_matrix<T> >
        >,
        result_of::scalar_product<Matrix, T>
    >::type
    prod (Matrix const & m, T const & t)
    {
        typedef typename result_of::scalar_product<Matrix, T>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, Matrix const &, T const &> ops;
        iterate<size<result_type> >(
            ops(retval, m, t), detail::matrix_scalar_mul_assign()
        );
        return retval;
    }

    template <typename Matrix, typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_matrix<Matrix>,
            mpl::not_<is_matrix<T> >
        >,
        result_of::scalar_product<Matrix, T>
    >::type
    prod (T const & t, Matrix const & m)
    { return prod(m, t); }

    template <typename Matrix, typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_matrix<Matrix>,
            mpl::not_<is_matrix<T> >
        >,
        result_of::scalar_quotient<Matrix, T>
    >::type
    div (Matrix const & m, T const & t)
    {
        typedef typename result_of::scalar_quotient<Matrix, T>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, Matrix const &, T const &> ops;
        iterate<size<result_type> >(
            ops(retval, m, t), detail::matrix_scalar_div_assign()
        );
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    template <typename MatrixL, typename MatrixR>
    typename lazy_enable_if<
        mpl::and_<
            is_matrix<MatrixL>,
            is_matrix<MatrixR>,
            mpl::equal_to<
                typename MatrixL::num_columns_t,
                typename MatrixR::num_rows_t
            >
        >,
        result_of::matrix_product<MatrixL, MatrixR>
    >::type
    operator* (MatrixL const & lhs, MatrixR const & rhs)
    { return prod(lhs, rhs); }

    template <typename Matrix, typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_matrix<Matrix>,
            mpl::not_<is_matrix<T> >
        >,
        result_of::scalar_product<Matrix, T>
    >::type
    operator* (Matrix const & m, T const & t)
    { return prod(m, t); }

    template <typename Matrix, typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_matrix<Matrix>,
            mpl::not_<is_matrix<T> >
        >,
        result_of::scalar_product<Matrix, T>
    >::type
    operator* (T const & t, Matrix const & m)
    { return prod(m, t); }

    template <typename Matrix, typename T>
    typename lazy_enable_if<
        mpl::and_<
            is_matrix<Matrix>,
            mpl::not_<is_matrix<T> >
        >,
        result_of::scalar_quotient<Matrix, T>
    >::type
    operator/ (Matrix const & m, T const & t)
    { return div(m, t); }

#endif

    template <typename MatrixL, typename MatrixR>
    typename lazy_enable_if<
        is_same_shape_matrix<MatrixL, MatrixR>,
        result_of::matrix_element_product<MatrixL, MatrixR>
    >::type
    element_prod (MatrixL const & lhs, MatrixR const & rhs)
    {
        typedef typename result_of::matrix_element_product<MatrixL, MatrixR>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, MatrixL const &, MatrixR const &> ops;
        iterate<size<result_type> >(
            ops(retval, lhs, rhs), detail::matrix_matrix_elem_mul_assign()
        );
        return retval;
    }

    template <typename MatrixL, typename MatrixR>
    typename lazy_enable_if<
        is_same_shape_matrix<MatrixL, MatrixR>,
        result_of::matrix_element_quotient<MatrixL, MatrixR>
    >::type
    element_div (MatrixL const & lhs, MatrixR const & rhs)
    {
        typedef typename result_of::matrix_element_quotient<MatrixL, MatrixR>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, MatrixL const &, MatrixR const &> ops;
        iterate<size<result_type> >(
            ops(retval, lhs, rhs), detail::matrix_matrix_elem_div_assign()
        );
        return retval;
    }

    template <typename Matrix>
    typename enable_if<
        is_matrix<Matrix>
    >::type
    swap (Matrix & lhs, Matrix & rhs)
    {
        typedef fusion::vector<Matrix &, Matrix &> ops;
        iterate<size<Matrix> >(
            ops(lhs, rhs), detail::swap()
        );
    }

    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        is_same_length_vector<VectorL, VectorR>,
        result_of::dot_product<VectorL, VectorR>
    >::type
    dot (VectorL const & lhs, VectorR const & rhs)
    {
        typedef typename result_of::detail::dot_product<
            typename fusion::result_of::value_at_c<typename result_of::transpose<VectorL>::type::value_types, 0>::type,
            typename fusion::result_of::value_at_c<typename result_of::transpose<VectorR>::type::value_types, 0>::type
        >::type result_type;
        result_type retval = detail::zero_value<result_type>::value();
        typedef fusion::vector<result_type &, VectorL const &, VectorR const &> ops;
        iterate<typename VectorL::num_rows_t>(
            ops(retval, lhs, rhs), detail::vector_vector_dot_product_assign()
        );
        return retval;
    }

    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        is_same_length_transpose_vector<VectorL, VectorR>,
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

    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        is_same_length_vector<VectorL, VectorR>,
        result_of::dot_product<VectorL, VectorR>
    >::type
    operator* (VectorL const & lhs, VectorR const & rhs)
    { return dot(lhs, rhs); }

    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        is_same_length_transpose_vector<VectorL, VectorR>,
        result_of::dot_product<VectorL, VectorR>
    >::type
    operator* (VectorL const & lhs, VectorR const & rhs)
    { return dot(lhs, rhs); }

#endif

    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_vector<VectorL, VectorR>,
            mpl::equal_to<
                typename VectorL::num_rows_t,
                size_t_<3>
            >
        >,
        result_of::cross_product<VectorL, VectorR>
    >::type
    cross (VectorL const & lhs, VectorR const & rhs)
    {
        typedef typename result_of::cross_product<VectorL, VectorR>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, VectorL const &, VectorR const &> ops;
        retval.at<0, 0>() = lhs.at<1, 0>() * rhs.at<2, 0>() - lhs.at<2, 0>() * rhs.at<1, 0>();
        retval.at<1, 0>() = lhs.at<2, 0>() * rhs.at<0, 0>() - lhs.at<0, 0>() * rhs.at<2, 0>();
        retval.at<2, 0>() = lhs.at<0, 0>() * rhs.at<1, 0>() - lhs.at<1, 0>() * rhs.at<0, 0>();
        return retval;
    }

    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_transpose_vector<VectorL, VectorR>,
            mpl::equal_to<
                typename VectorL::num_columns_t,
                size_t_<3>
            >
        >,
        result_of::cross_product<VectorL, VectorR>
    >::type
    cross (VectorL const & lhs, VectorR const & rhs)
    {
        typedef typename result_of::cross_product<VectorL, VectorR>::type result_type;
        result_type retval;
        typedef fusion::vector<result_type &, VectorL const &, VectorR const &> ops;
        retval.at<0, 0>() = lhs.at<0, 1>() * rhs.at<0, 2>() - lhs.at<0, 2>() * rhs.at<0, 1>();
        retval.at<0, 1>() = lhs.at<0, 2>() * rhs.at<0, 0>() - lhs.at<0, 0>() * rhs.at<0, 2>();
        retval.at<0, 2>() = lhs.at<0, 0>() * rhs.at<0, 1>() - lhs.at<0, 1>() * rhs.at<0, 0>();
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_vector<VectorL, VectorR>,
            mpl::equal_to<
                typename VectorL::num_rows_t,
                size_t_<3>
            >
        >,
        result_of::cross_product<VectorL, VectorR>
    >::type
    operator^ (VectorL const & lhs, VectorR const & rhs)
    { return cross(lhs, rhs); }

    template <typename VectorL, typename VectorR>
    typename lazy_enable_if<
        mpl::and_<
            is_same_length_transpose_vector<VectorL, VectorR>,
            mpl::equal_to<
                typename VectorL::num_columns_t,
                size_t_<3>
            >
        >,
        result_of::cross_product<VectorL, VectorR>
    >::type
    operator^ (VectorL const & lhs, VectorR const & rhs)
    { return cross(lhs, rhs); }

#endif

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

    template <typename Vector>
    typename lazy_enable_if<
        is_transpose_vector<Vector>,
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

    template <typename Vector>
    typename lazy_enable_if<
        is_transpose_vector<Vector>,
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

    template <typename Vector>
    typename lazy_enable_if<
        is_transpose_vector<Vector>,
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

    template <typename Vector>
    typename lazy_enable_if<
        is_vector<Vector>,
        result_of::sum<Vector>
    >::type
    norm_inf (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        typedef std::pair<result_type, std::size_t> temp_type;
        temp_type tmp(detail::abs_(static_cast<result_type>(v.at<0, 0>())), 0);
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_rows_t>(
            ops(tmp, v), detail::vector_norm_inf()
        );
        return tmp.first;
    }

    template <typename Vector>
    typename lazy_enable_if<
        is_transpose_vector<Vector>,
        result_of::sum<Vector>
    >::type
    norm_inf (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type result_type;
        typedef std::pair<result_type, std::size_t> temp_type;
        temp_type tmp(detail::abs_(static_cast<result_type>(v.at<0, 0>())), 0);
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_columns_t>(
            ops(tmp, v), detail::transpose_vector_norm_inf()
        );
        return tmp.first;
    }

    template <typename Vector>
    typename enable_if<
        is_vector<Vector>,
        std::size_t
    >::type
    norm_inf_index (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type sum_type;
        typedef std::pair<sum_type, std::size_t> temp_type;
        temp_type tmp(detail::abs_(static_cast<sum_type>(v.at<0, 0>())), 0);
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_rows_t>(
            ops(tmp, v), detail::vector_norm_inf()
        );
        return tmp.second;
    }

    template <typename Vector>
    typename enable_if<
        is_transpose_vector<Vector>,
        std::size_t
    >::type
    norm_inf_index (Vector const & v)
    {
        typedef typename result_of::sum<Vector>::type sum_type;
        typedef std::pair<sum_type, std::size_t> temp_type;
        temp_type tmp(detail::abs_(static_cast<sum_type>(v.at<0, 0>())), 0);
        typedef fusion::vector<temp_type &, Vector const &> ops;
        iterate<typename Vector::num_columns_t>(
            ops(tmp, v), detail::transpose_vector_norm_inf()
        );
        return tmp.second;
    }

    // TODO: Implement these for expression template optimizations, as needed.
    // axpy_prod(A, u, w, true);  // w = A * u
    // axpy_prod(A, u, w, false); // w += A * u
    // axpy_prod(u, A, w, true);  // w = trans(A) * u
    // axpy_prod(u, A, w, false); // w += trans(A) * u
    // axpy_prod(A, B, C, true);  // C = A * B
    // axpy_prod(A, B, C, false); // C += A * B

    // TODO: When are these preferable to the above?
    // opb_prod(A, B, C, true);  // C = A * B
    // opb_prod(A, B, C, false); // C += A * B

    template <typename Matrix>
    typename lazy_enable_if<
        mpl::and_<
            is_square_matrix<Matrix>,
            mpl::equal_to<
                typename Matrix::num_rows_t,
                size_t_<1>
            >
        >,
        result_of::determinant<Matrix>
    >::type
    determinant (Matrix const & m)
    { return m.at<0, 0>(); }

    template <typename Matrix>
    typename lazy_enable_if<
        mpl::and_<
            is_square_matrix<Matrix>,
            mpl::equal_to<
                typename Matrix::num_rows_t,
                size_t_<2>
            >
        >,
        result_of::determinant<Matrix>
    >::type
    determinant (Matrix const & m)
    { return m.at<0, 0>() * m.at<1, 1>() - m.at<0, 1>() * m.at<1, 0>(); }

    template <typename Matrix>
    typename lazy_enable_if<
        mpl::and_<
            is_square_matrix<Matrix>,
            mpl::equal_to<
                typename Matrix::num_rows_t,
                size_t_<3>
            >
        >,
        result_of::determinant<Matrix>
    >::type
    determinant (Matrix const & m)
    {
        return
            m.at<0, 0>() * (m.at<1, 1>() * m.at<2, 2>() - m.at<1, 2>() * m.at<2, 1>()) -
            m.at<0, 1>() * (m.at<1, 0>() * m.at<2, 2>() - m.at<1, 2>() * m.at<2, 0>()) +
            m.at<0, 2>() * (m.at<1, 0>() * m.at<2, 1>() - m.at<1, 1>() * m.at<2, 0>());
    }

    template <typename Matrix>
    typename lazy_enable_if<
        mpl::and_<
            is_square_matrix<Matrix>,
            mpl::less<
                size_t_<3>,
                typename Matrix::num_rows_t
            >
        >,
        result_of::determinant<Matrix>
    >::type
    determinant (Matrix const & m)
    {
        typedef typename result_of::determinant<Matrix>::type result_type;
        result_type retval = detail::one_value<result_type>::value();
        typedef typename detail::get_value_type<result_type>::type raw_value_type;
        typedef array<
            array<raw_value_type, Matrix::num_columns_t::value>,
            Matrix::num_rows_t::value
        > temp_matrix_type;
        temp_matrix_type temp_matrix;
        typedef fusion::vector<temp_matrix_type &, Matrix const &> ops1;
        iterate<size<Matrix> >(
            ops1(temp_matrix, m), detail::matrix_to_temp_assign<result_type>()
        );
        array<std::size_t, Matrix::num_rows_t::value> indices;
        try {
            retval *= detail::lu_decompose(temp_matrix, indices);
            typedef fusion::vector<result_type &, temp_matrix_type const &> ops2;
            iterate<typename Matrix::num_rows_t>(
                ops2(retval, temp_matrix), detail::accumulate_determinant()
            );
        } catch (singular_matrix const &) {
            retval = detail::zero_value<result_type>::value();
        }
        return retval;
    }

    template <typename Matrix>
    typename lazy_enable_if<
        is_square_matrix<Matrix>,
        result_of::inverse<Matrix>
    >::type
    inverse (Matrix const & m)
    {
        typedef typename result_of::inverse<Matrix>::type result_type;
        result_type retval;
        typedef typename result_of::determinant<Matrix>::type temp_value_type;
        typedef typename detail::get_value_type<temp_value_type>::type raw_value_type;
        typedef array<
            array<raw_value_type, Matrix::num_columns_t::value>,
            Matrix::num_rows_t::value
        > temp_matrix_type;
        temp_matrix_type temp_matrix;
        typedef fusion::vector<temp_matrix_type &, Matrix const &> ops1;
        iterate<size<Matrix> >(
            ops1(temp_matrix, m), detail::matrix_to_temp_assign<temp_value_type>()
        );
        typedef array<std::size_t, Matrix::num_rows_t::value> indices_type;
        indices_type indices;
        detail::lu_decompose(temp_matrix, indices);
        temp_matrix_type temp_result_matrix;
        typedef fusion::vector<temp_matrix_type &, temp_matrix_type const &, indices_type const &> ops2;
        iterate<typename Matrix::num_rows_t>(
            ops2(temp_result_matrix, temp_matrix, indices), detail::assign_inverted_column()
        );
        typedef fusion::vector<result_type &, temp_matrix_type const &> ops3;
        iterate<size<Matrix> >(
            ops3(retval, temp_result_matrix), detail::temp_to_matrix_assign()
        );
        return retval;
    }

    template <typename AMatrix, typename XVector, typename BVector>
    typename enable_if<
        mpl::and_<
            is_square_matrix<AMatrix>,
            is_same_length_vector<XVector, BVector>,
            mpl::equal_to<
                typename AMatrix::num_columns_t,
                typename XVector::num_rows_t
            >
        >
    >::type
    solve (AMatrix const & A, BVector const & b, XVector & x)
    {
        typedef BOOST_TYPEOF((AMatrix() * XVector())) a_times_b_type;
        BOOST_MPL_ASSERT((is_convertible<a_times_b_type, BVector>));

        typedef typename result_of::determinant<AMatrix>::type temp_value_type;
        typedef typename detail::get_value_type<temp_value_type>::type raw_value_type;
        typedef array<
            array<raw_value_type, AMatrix::num_columns_t::value>,
            AMatrix::num_rows_t::value
        > temp_matrix_type;
        temp_matrix_type temp_A;
        typedef fusion::vector<temp_matrix_type &, AMatrix const &> ops1;
        iterate<size<AMatrix> >(
            ops1(temp_A, A), detail::matrix_to_temp_assign<temp_value_type>()
        );
        array<std::size_t, AMatrix::num_rows_t::value> indices;
        detail::lu_decompose(temp_A, indices);
        typedef array<raw_value_type, AMatrix::num_rows_t::value> temp_vector_type;
        temp_vector_type temp_vector;
        typedef fusion::vector<temp_vector_type &, BVector const &> ops2;
        iterate<typename AMatrix::num_rows_t>(
            ops2(temp_vector, b), detail::matrix_to_temp_vector_assign<temp_value_type>()
        );
        detail::lu_substitute(temp_A, indices, temp_vector);
        typedef fusion::vector<XVector &, temp_vector_type const &> ops3;
        iterate<typename XVector::num_rows_t>(
            ops3(x, temp_vector), detail::temp_vector_to_matrix_assign()
        );
    }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_OPERATIONS_HPP
