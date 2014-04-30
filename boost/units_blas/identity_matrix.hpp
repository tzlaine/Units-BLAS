// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_IDENTIY_MATRIX_HPP
#define BOOST_UNITS_BLAS_IDENTIY_MATRIX_HPP

#include <boost/units_blas/matrix.hpp>
#include <boost/units_blas/detail/inverse_type.hpp>
#include <boost/units_blas/detail/one_value.hpp>
#include <boost/units_blas/detail/zero_value.hpp>


namespace boost { namespace units_blas {

    namespace detail {

        template <typename Matrix>
        struct identity_assign
        {
            template <std::size_t I>
            void call ()
            {
                constexpr std::size_t row = I / Matrix::num_columns;
                constexpr std::size_t column = I % Matrix::num_columns;
                using type = tuple_element_t<I, Matrix>;
                tuple_access::get<I>(m_) =
                    row == column ?
                    one_value<type>::value() :
                    zero_value<type>::value();
            }

            Matrix & m_;
        };

    }

    /** Returns the type that can be used as the left-identity matrix for @c
        Matrix.  Note that not every matrix type has an identiy matrix
        type. */
    template <typename Matrix>
    struct left_identity_matrix_type
    {
        using type = decltype(
            prod(
                std::declval<Matrix>(),
                std::declval<detail::inverse_type_t<Matrix>>()
            )
        );
    };

    /** Returns the type that can be used as the right-identity matrix for @c
        Matrix.  Note that not every matrix type has an identiy matrix
        type. */
    template <typename Matrix>
    struct right_identity_matrix_type
    {
        using type = decltype(
            prod(
                std::declval<detail::inverse_type_t<Matrix>>(),
                std::declval<Matrix>()
            )
        );
    };

    template <typename Matrix>
    using left_identity_matrix =
        typename left_identity_matrix_type<Matrix>::type;

    template <typename Matrix>
    using right_identity_matrix =
        typename right_identity_matrix_type<Matrix>::type;

    /** Returns a matrix<> I whose diagonal elements are 1, and whose
        nondiagonal elements are 0, and whose element types are such that for
        any @c Matrix m, type of the expression m * I is assignable to a @c
        Matrix.  Further, for all i,j in [0, @c Matrix::num_rows_t::value), (m
        * I).at<i, j>() is within epsilon of m.at<i, j>().  Note that not
        every matrix type has an left-identity matrix type. */
    template <typename Matrix>
    auto make_left_identity_matrix ()
    {
        using result_type = left_identity_matrix<Matrix>;
        result_type retval;
        detail::iterate_simple<result_type::num_elements>(
            detail::identity_assign<result_type>{retval}
        );
        return retval;
    }

    /** Returns a matrix<> I whose diagonal elements are 1, and whose
        nondiagonal elements are 0, and whose element types are such that for
        any @c Matrix m, type of the expression m * I is assignable to a @c
        Matrix.  Further, for all i,j in [0, @c Matrix::num_rows_t::value), (m
        * I).at<i, j>() is within epsilon of m.at<i, j>().  Note that not
        every matrix type has an right-identity matrix type. */
    template <typename Matrix>
    auto make_right_identity_matrix ()
    {
        using result_type = right_identity_matrix<Matrix>;
        result_type retval;
        detail::iterate_simple<result_type::num_elements>(
            detail::identity_assign<result_type>{retval}
        );
        return retval;
    }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MATRIX_HPP
