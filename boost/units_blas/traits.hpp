// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_TRAITS_HPP
#define BOOST_UNITS_BLAS_TRAITS_HPP

#include <boost/units_blas/matrix_fwd.hpp>


namespace boost { namespace units_blas {

    /** Metafunction that returns the number of elements in \a Matrix. */
    template <typename Matrix>
    struct size :
        mpl::times<typename Matrix::num_rows_t, typename Matrix::num_columns_t>::type
    {};

    /** Metafunction that returns the number of rows in \a Matrix. */
    template <typename Matrix>
    struct rows : Matrix::num_rows_t {};

    /** Metafunction that returns the number of columns in \a Matrix. */
    template <typename Matrix>
    struct columns : Matrix::num_columns_t {};

    /** Metafunction that returns true iff \a T is a (possibly cv-qualified)
        matrix<>. */
    template <typename T>
    struct is_matrix : mpl::false_ {};
#ifndef BOOST_UNITS_BLAS_DOXYGEN
    template <typename T>
    struct is_matrix<matrix<T> > : mpl::true_ {};
    template <typename T>
    struct is_matrix<matrix<T> const> : mpl::true_ {};
    template <typename T>
    struct is_matrix<matrix<T> volatile> : mpl::true_ {};
    template <typename T>
    struct is_matrix<matrix<T> const volatile> : mpl::true_ {};
#endif

    /** Metafunction that returns true iff \a T is a (possibly cv-qualified)
        matrix<> with the same number of rows as columns. */
    template <typename T>
    struct is_square_matrix :
        mpl::and_<
            is_matrix<T>,
            mpl::equal_to<rows<T>, columns<T> >
        >::type
    {};

    /** Metafunction that returns true iff \a T is a (possibly cv-qualified)
        matrix<> with exactly 1 column. */
    template <typename T>
    struct is_vector :
        mpl::and_<
            is_matrix<T>,
            mpl::equal_to<columns<T>, size_t_<1> >
        >::type
    {};

    /** Metafunction that returns true iff \a T is a (possibly cv-qualified)
        matrix<> with exactly 1 row. */
    template <typename T>
    struct is_transpose_vector :
        mpl::and_<
            is_matrix<T>,
            mpl::equal_to<rows<T>, size_t_<1> >
        >::type
    {};

    /** Metafunction that returns true iff \a T0 and \a T1 are both (possibly
        cv-qualified) matrix<>s with the same dimensions. */
    template <typename T0, typename T1>
    struct is_same_shape_matrix :
        mpl::and_<
            is_matrix<T0>,
            is_matrix<T1>,
            mpl::equal_to<rows<T0>, rows<T1> >,
            mpl::equal_to<columns<T0>, columns<T1> >
        >::type
    {};

    /** Metafunction that returns true iff \a T0 and \a T1 are both (possibly
        cv-qualified) "vector" matrix<>s with the same dimensions. */
    template <typename T0, typename T1>
    struct is_same_length_vector :
        mpl::and_<
            is_vector<T0>,
            is_vector<T1>,
            mpl::equal_to<rows<T0>, rows<T1> >
        >::type
    {};

    /** Metafunction that returns true iff \a T0 and \a T1 are both (possibly
        cv-qualified) "transpose vector" matrix<>s with the same dimensions. */
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
