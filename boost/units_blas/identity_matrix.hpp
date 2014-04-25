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
                using type = typename std::tuple_element<
                    I,
                    typename Matrix::value_types
                >::type;
                tuple_access::get<I>(m_) =
                    row == column ?
                    one_value<type>::value() :
                    zero_value<type>::value();
            }

            Matrix & m_;
        };

    }

    /** Returns the type that can be used as the identity matrix for @c
        Matrix.  Note that not every matrix type has an identiy matrix
        type. */
    template <typename Matrix>
    struct identity_matrix
    {
        using type = decltype(
            std::declval<Matrix>() *
            std::declval<typename detail::inverse_type<Matrix>::type>()
        );
    };

    /** Returns a matrix<> I whose diagonal elements are 1, and whose
        nondiagonal elements are 0, and whose element types are such that for
        any @c Matrix m, type of the expression m * I is assignable to a @c
        Matrix.  Further, for all i,j in [0, @c Matrix::num_rows_t::value), (m
        * I).at<i, j>() is within epsilon of m.at<i, j>().  Note that not
        every matrix type has an identiy matrix type. */
    template <typename Matrix>
    auto make_identity_matrix ()
    {
        using result_type = typename identity_matrix<Matrix>::type;
        result_type retval;
        detail::iterate_simple<result_type::num_elements>(
            detail::identity_assign<result_type>{retval}
        );
        return retval;
    }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MATRIX_HPP
