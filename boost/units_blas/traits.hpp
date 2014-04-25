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
#include <type_traits>


namespace boost { namespace units_blas {

    template <typename T>
    struct is_matrix :
        std::false_type
    {};

#ifndef BOOST_UNITS_BLAS_DOXYGEN
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    struct is_matrix<matrix_t<Tuple, Rows, Columns>> :
        std::true_type
    {};
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    struct is_matrix<matrix_t<Tuple, Rows, Columns> const> :
        std::true_type
    {};
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    struct is_matrix<matrix_t<Tuple, Rows, Columns> volatile> :
        std::true_type
    {};
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    struct is_matrix<matrix_t<Tuple, Rows, Columns> const volatile> :
        std::true_type
    {};

    namespace detail {

        std::false_type is_matrix_or_derived_impl(...);

        template <typename Tuple, std::size_t Rows, std::size_t Columns>
        std::true_type is_matrix_or_derived_impl(
            matrix_t<Tuple, Rows, Columns>
        );

    }
#endif

    template <typename T>
    using is_matrix_t = typename is_matrix<T>::type;

#if 0 // TODO: Enable when Clang fixes the variable template internal linkage bug.
    template <typename T>
    const bool is_matrix_v = is_matrix<T>::value;
#endif

    /** Returns true iff cv-unqualified @c T is, or is derived from, a @c
        matrix<>. */
    template <typename T>
    struct is_matrix_or_derived :
        decltype(detail::is_matrix_or_derived_impl(
            std::declval<std::remove_cv_t<T>>()
        ))
    {};

    template <typename T>
    using is_matrix_or_derived_t = typename is_matrix_or_derived<T>::type;

#if 0 // TODO: Enable when Clang fixes the variable template internal linkage bug.
    template <typename T>
    const bool is_matrix_or_derived_v = is_matrix_or_derived<T>::value;
#endif

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_TRAITS_HPP
