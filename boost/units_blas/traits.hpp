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

    /** Returns the number of elements in @c Matrix. */
    template <typename Matrix>
    struct size :
        mpl::times<typename Matrix::num_rows_t, typename Matrix::num_columns_t>::type
    {};

    /** Returns the number of rows in @c Matrix. */
    template <typename Matrix>
    struct rows : Matrix::num_rows_t {};

    /** Returns the number of columns in @c Matrix. */
    template <typename Matrix>
    struct columns : Matrix::num_columns_t {};

    /** Returns true iff @c T is a (possibly cv-qualified) @c matrix<>. */
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

    /** Returns true iff @c T is a (possibly cv-qualified) @c matrix<> with
        the same number of rows as columns. */
    template <typename T>
    struct is_square_matrix :
        mpl::and_<
            is_matrix<T>,
            mpl::equal_to<rows<T>, columns<T> >
        >::type
    {};

    /** Returns true iff @c T is a (possibly cv-qualified) @c matrix<> with
        exactly 1 column. */
    template <typename T>
    struct is_vector :
        mpl::and_<
            is_matrix<T>,
            mpl::equal_to<columns<T>, mpl::size_t<1> >
        >::type
    {};

    /** Returns true iff @c T is a (possibly cv-qualified) @c matrix<> with
        exactly 1 row. */
    template <typename T>
    struct is_transpose_vector :
        mpl::and_<
            is_matrix<T>,
            mpl::equal_to<rows<T>, mpl::size_t<1> >
        >::type
    {};

    /** Returns true iff @c T0 and @c T1 are both (possibly cv-qualified)
        <c>matrix<></c>s with the same dimensions. */
    template <typename T0, typename T1>
    struct is_same_shape_matrix :
        mpl::and_<
            is_matrix<T0>,
            is_matrix<T1>,
            mpl::equal_to<rows<T0>, rows<T1> >,
            mpl::equal_to<columns<T0>, columns<T1> >
        >::type
    {};

    /** Returns true iff @c T0 and @c T1 are both (possibly cv-qualified)
        "vector" <c>matrix<></c>s with the same dimensions. */
    template <typename T0, typename T1>
    struct is_same_length_vector :
        mpl::and_<
            is_vector<T0>,
            is_vector<T1>,
            mpl::equal_to<rows<T0>, rows<T1> >
        >::type
    {};

    /** Returns true iff @c T0 and @c T1 are both (possibly cv-qualified)
        "transpose vector" <c>matrix<></c>s with the same dimensions. */
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
