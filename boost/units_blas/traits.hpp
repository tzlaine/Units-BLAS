// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_TRAITS_HPP
#define BOOST_UNITS_BLAS_TRAITS_HPP

#include <boost/units_blas/matrix.hpp>


namespace boost { namespace units_blas {

    template <typename Matrix>
    struct size :
        mpl::times<typename Matrix::num_rows_t, typename Matrix::num_columns_t>::type
    {};

    template <typename Matrix>
    struct rows : Matrix::num_rows_t {};

    template <typename Matrix>
    struct columns : Matrix::num_columns_t {};

    template <typename T>
    struct is_matrix : mpl::false_ {};
    template <typename T>
    struct is_matrix<matrix<T> > : mpl::true_ {};
    template <typename T>
    struct is_matrix<matrix<T> const> : mpl::true_ {};
    template <typename T>
    struct is_matrix<matrix<T> volatile> : mpl::true_ {};
    template <typename T>
    struct is_matrix<matrix<T> const volatile> : mpl::true_ {};

    template <typename T>
    struct is_square_matrix :
        mpl::and_<
            is_matrix<T>,
            mpl::equal_to<rows<T>, columns<T> >
        >::type
    {};

    template <typename T>
    struct is_vector :
        mpl::and_<
            is_matrix<T>,
            mpl::equal_to<columns<T>, size_t_<1> >
        >::type
    {};

    template <typename T>
    struct is_transpose_vector :
        mpl::and_<
            is_matrix<T>,
            mpl::equal_to<rows<T>, size_t_<1> >
        >::type
    {};

    template <typename T0, typename T1>
    struct is_same_shape_matrix :
        mpl::and_<
            is_matrix<T0>,
            is_matrix<T1>,
            mpl::equal_to<rows<T0>, rows<T1> >,
            mpl::equal_to<columns<T0>, columns<T1> >
        >::type
    {};

    template <typename T0, typename T1>
    struct is_same_length_vector :
        mpl::and_<
            is_vector<T0>,
            is_vector<T1>,
            mpl::equal_to<rows<T0>, rows<T1> >
        >::type
    {};

    template <typename T0, typename T1>
    struct is_same_length_transpose_vector :
        mpl::and_<
            is_transpose_vector<T0>,
            is_transpose_vector<T1>,
            mpl::equal_to<columns<T0>, columns<T1> >
        >::type
    {};

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_TRAITS_HPP
