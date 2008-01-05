// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DETAIL_SIMPLE_ITERATION_HPP
#define BOOST_UNITS_BLAS_DETAIL_SIMPLE_ITERATION_HPP

#include <boost/units_blas/iterate.hpp>
#include <boost/units_blas/detail/lu.hpp>
#include <boost/units_blas/detail/one_value.hpp>

#include <boost/fusion/sequence/intrinsic/at.hpp>


namespace boost { namespace units_blas { namespace detail {

    struct assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                std::size_t const I = N / to_matrix::num_columns_t::value;
                std::size_t const J = N % to_matrix::num_columns_t::value;
                typedef typename result_of::value_at_c<to_matrix, I, J>::type to_type;
                fusion::at_c<0>(operands).template at<I, J>() =
                    static_cast<to_type>(fusion::at_c<1>(operands).template at<I, J>());
            }
    };

    struct plus_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                std::size_t const I = N / to_matrix::num_columns_t::value;
                std::size_t const J = N % to_matrix::num_columns_t::value;
                typedef typename result_of::value_at_c<to_matrix, I, J>::type to_type;
                fusion::at_c<0>(operands).template at<I, J>() =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands).template at<I, J>() +
                        fusion::at_c<1>(operands).template at<I, J>()
                    );
            }
    };

    struct minus_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                std::size_t const I = N / to_matrix::num_columns_t::value;
                std::size_t const J = N % to_matrix::num_columns_t::value;
                typedef typename result_of::value_at_c<to_matrix, I, J>::type to_type;
                fusion::at_c<0>(operands).template at<I, J>() =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands).template at<I, J>() -
                        fusion::at_c<1>(operands).template at<I, J>()
                    );
            }
    };

    struct matrix_scalar_mul_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                std::size_t const I = N / to_matrix::num_columns_t::value;
                std::size_t const J = N % to_matrix::num_columns_t::value;
                typedef typename result_of::value_at_c<to_matrix, I, J>::type to_type;
                fusion::at_c<0>(operands).template at<I, J>() =
                    static_cast<to_type>(
                        fusion::at_c<1>(operands).template at<I, J>() *
                        fusion::at_c<2>(operands)
                    );
            }
    };

    struct matrix_scalar_div_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                std::size_t const I = N / to_matrix::num_columns_t::value;
                std::size_t const J = N % to_matrix::num_columns_t::value;
                typedef typename result_of::value_at_c<to_matrix, I, J>::type to_type;
                fusion::at_c<0>(operands).template at<I, J>() =
                    static_cast<to_type>(
                        fusion::at_c<1>(operands).template at<I, J>() /
                        fusion::at_c<2>(operands)
                    );
            }
    };

    struct identity_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                typedef typename result_of::value_at_c<to_matrix, N, N>::type to_type;
                fusion::at_c<0>(operands).template at<N, N>() = one_value<to_type>::value();
            }
    };

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_DETAIL_SIMPLE_ITERATION_HPP
