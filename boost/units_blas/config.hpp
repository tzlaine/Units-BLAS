// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_CONFIG_HPP
#define BOOST_UNITS_BLAS_CONFIG_HPP


#ifndef BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS
    /** This is an overridable configuration macro controlling whether operator
        overloads are used in units_blas.  See the detailed Configuration page
        of the documentation for details.  To disable the use of operator
        overloads, define this macro to be zero before including any library
        header.*/
    #define BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS 1
#endif

#ifndef BOOST_UNITS_BLAS_USE_INEXACT_DETERMINANT_TYPE
    /** This is an overridable configuration macro controlling whether fully
        accurate, but slow (O(n!) in the number of rows), or inexact but fast
        (O(n) in the number of rows) metafunction should be used to determine
        the type of a determinant.  See the detailed Configuration page of the
        documentation for details.  To use the inexact version, define this
        macro to a nonzero value before including any library header.*/
    #define BOOST_UNITS_BLAS_USE_INEXACT_DETERMINANT_TYPE 0
#endif

#ifndef BOOST_UNITS_BLAS_DOXYGEN

#define BOOST_UNITS_BLAS_INLINE inline
#ifdef NDEBUG
    #if defined(_MSC_VER)
        #undef BOOST_UNITS_BLAS_INLINE
        #define BOOST_UNITS_BLAS_INLINE __forceinline
    #elif defined(__GNUC__) && __GNUC__ > 3
        #undef BOOST_UNITS_BLAS_INLINE
        #define BOOST_UNITS_BLAS_INLINE inline __attribute__ ((always_inline))
    #endif
#endif

#endif

#ifdef BOOST_UNITS_BLAS_DOXYGEN
namespace detail {
    struct unspecified {};
}
#endif

#endif // BOOST_UNITS_BLAS_CONFIG_HPP
