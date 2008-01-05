// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DETAIL_ITERATION_HPP
#define BOOST_UNITS_BLAS_DETAIL_ITERATION_HPP

#include <boost/units_blas/iterate.hpp>
#include <boost/units_blas/detail/lu.hpp>

#include <boost/units/quantity.hpp>


namespace boost { namespace units_blas { namespace detail {

    struct negate_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type matrix;
                std::size_t const I = N / matrix::num_columns_t::value;
                std::size_t const J = N % matrix::num_columns_t::value;
                fusion::at_c<0>(operands).template at<I, J>() =
                    -fusion::at_c<1>(operands).template at<I, J>();
            }
    };

    template <typename Rows, typename Columns>
    struct slice_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                std::size_t const I = N / to_matrix::num_columns_t::value;
                std::size_t const J = N % to_matrix::num_columns_t::value;
                fusion::at_c<0>(operands).template at<I, J>() =
                    fusion::at_c<1>(operands).template at<
                        mpl::at_c<Rows, I>::type::value,
                        mpl::at_c<Columns, J>::type::value
                    >();
            }
    };

    struct transpose_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                std::size_t const I = N / to_matrix::num_columns_t::value;
                std::size_t const J = N % to_matrix::num_columns_t::value;
                typedef typename result_of::value_at_c<to_matrix, I, J>::type to_type;
                fusion::at_c<0>(operands).template at<I, J>() =
                    static_cast<to_type>(fusion::at_c<1>(operands).template at<J, I>());
            }
    };

    struct matrix_matrix_mul_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                typedef typename operand_n<Operands, 1>::type from_matrix_lhs;
                std::size_t const I = N / (from_matrix_lhs::num_columns_t::value * to_matrix::num_columns_t::value);
                std::size_t const J = (N / from_matrix_lhs::num_columns_t::value) % to_matrix::num_columns_t::value;
                std::size_t const K = N % from_matrix_lhs::num_columns_t::value;
                typedef typename result_of::value_at_c<to_matrix, I, J>::type to_type;
                fusion::at_c<0>(operands).template at<I, J>() =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands).template at<I, J>() +
                        fusion::at_c<1>(operands).template at<I, K>() *
                        fusion::at_c<2>(operands).template at<K, J>()
                    );
            }
    };

    struct matrix_matrix_elem_add_assign
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
                        fusion::at_c<1>(operands).template at<I, J>() +
                        fusion::at_c<2>(operands).template at<I, J>()
                    );
            }
    };

    struct matrix_matrix_elem_sub_assign
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
                        fusion::at_c<1>(operands).template at<I, J>() -
                        fusion::at_c<2>(operands).template at<I, J>()
                    );
            }
    };

    struct matrix_matrix_elem_mul_assign
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
                        fusion::at_c<2>(operands).template at<I, J>()
                    );
            }
    };

    struct matrix_matrix_elem_div_assign
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
                        fusion::at_c<2>(operands).template at<I, J>()
                    );
            }
    };

    struct swap
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type matrix;
                std::size_t const I = N / matrix::num_columns_t::value;
                std::size_t const J = N % matrix::num_columns_t::value;
                typedef typename result_of::value_at_c<matrix, I, J>::type value_type;
                value_type tmp = fusion::at_c<0>(operands).template at<I, J>();
                fusion::at_c<0>(operands).template at<I, J>() =
                    fusion::at_c<1>(operands).template at<I, J>();
                fusion::at_c<1>(operands).template at<I, J>() = tmp;
            }
    };

    struct vector_vector_dot_product_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_type;
                fusion::at_c<0>(operands) =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands) +
                        fusion::at_c<1>(operands).template at<N, 0>() *
                        fusion::at_c<2>(operands).template at<N, 0>()
                    );
            }
    };

    struct transpose_vector_transpose_vector_dot_product_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_type;
                fusion::at_c<0>(operands) =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands) +
                        fusion::at_c<1>(operands).template at<0, N>() *
                        fusion::at_c<2>(operands).template at<0, N>()
                    );
            }
    };

    template <bool absolute_value>
    struct vector_sum
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_type;
                const typename result_of::value_at_c<
                    typename operand_n<Operands, 1>::type, N, 0
                >::type &val = fusion::at_c<1>(operands).template at<N, 0>();
                fusion::at_c<0>(operands) =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands) +
                        (absolute_value ? abs_(val) : val)
                    );
            }
    };

    template <bool absolute_value>
    struct transpose_vector_sum
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_type;
                const typename result_of::value_at_c<
                    typename operand_n<Operands, 1>::type, 0, N
                >::type &val = fusion::at_c<1>(operands).template at<0, N>();
                fusion::at_c<0>(operands) =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands) +
                        (absolute_value ? abs_(val) : val)
                    );
            }
    };

    struct vector_norm_2
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_type;
                fusion::at_c<0>(operands) =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands) +
                        fusion::at_c<1>(operands).template at<N, 0>() *
                        fusion::at_c<1>(operands).template at<N, 0>()
                    );
            }
    };

    struct transpose_vector_norm_2
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_type;
                fusion::at_c<0>(operands) =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands) +
                        fusion::at_c<1>(operands).template at<0, N>() *
                        fusion::at_c<1>(operands).template at<0, N>()
                    );
            }
    };

    struct vector_norm_inf
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type::first_type compare_type;
                compare_type val =
                    abs_(static_cast<compare_type>(fusion::at_c<1>(operands).template at<N, 0>()));
                if (fusion::at_c<0>(operands).first <= val) {
                    fusion::at_c<0>(operands).first = val;
                    fusion::at_c<0>(operands).second = N;
                }
            }
    };

    struct transpose_vector_norm_inf
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type::first_type compare_type;
                compare_type val =
                    abs_(static_cast<compare_type>(fusion::at_c<1>(operands).template at<0, N>()));
                if (fusion::at_c<0>(operands).first <= val) {
                    fusion::at_c<0>(operands).first = val;
                    fusion::at_c<0>(operands).second = N;
                }
            }
    };

    template <typename T>
    struct matrix_to_temp_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 1>::type from_matrix;
                std::size_t const I = N / from_matrix::num_columns_t::value;
                std::size_t const J = N % from_matrix::num_columns_t::value;
                typedef typename operand_n<Operands, 0>::type::value_type::value_type to_type;
                fusion::at_c<0>(operands)[I][J] =
                    static_cast<to_type>(fusion::at_c<1>(operands).template at<I, J>());
            }
    };

    template <typename Unit, typename ValueType>
    struct matrix_to_temp_assign<units::quantity<Unit, ValueType> >
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 1>::type from_matrix;
                std::size_t const I = N / from_matrix::num_columns_t::value;
                std::size_t const J = N % from_matrix::num_columns_t::value;
                typedef typename operand_n<Operands, 0>::type::value_type::value_type to_type;
                fusion::at_c<0>(operands)[I][J] =
                    static_cast<to_type>(fusion::at_c<1>(operands).template at<I, J>().value());
            }
    };

    template <typename T>
    struct assign_from_value_type
    {
        template <typename ValueType>
        static void call (T & x, ValueType const & value)
            { x = static_cast<T>(value); }
    };

    template <typename Unit, typename ValueType>
    struct assign_from_value_type<units::quantity<Unit, ValueType> >
    {
        template <typename ValueType2>
        static void call (units::quantity<Unit, ValueType> & x, ValueType2 const & value)
            { x = units::quantity<Unit, ValueType>::from_value(static_cast<ValueType>(value)); }
    };

    struct temp_to_matrix_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                std::size_t const I = N / to_matrix::num_columns_t::value;
                std::size_t const J = N % to_matrix::num_columns_t::value;
                typedef typename result_of::value_at_c<to_matrix, I, J>::type to_type;
                assign_from_value_type<to_type>::call(
                    fusion::at_c<0>(operands).template at<I, J>(),
                    fusion::at_c<1>(operands)[I][J]
                );
            }
    };

    struct accumulate_determinant
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_type;
                fusion::at_c<0>(operands) =
                    static_cast<to_type>(
                        fusion::at_c<0>(operands) *
                        fusion::at_c<1>(operands)[N][N]
                    );
            }
    };

    struct assign_inverted_column
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type::value_type::value_type value_type;
                std::size_t const size = operand_n<Operands, 0>::type::static_size;
                array<value_type, size> columns;
                columns.assign(value_type(0));
                columns[N] = value_type(1.0);
                detail::lu_substitute(fusion::at_c<1>(operands), fusion::at_c<2>(operands), columns);
                for (std::size_t i = 0; i < size; ++i) {
                    fusion::at_c<0>(operands)[i][N] = columns[i];
                }
            }
    };

    template <typename T>
    struct matrix_to_temp_vector_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type::value_type to_type;
                fusion::at_c<0>(operands)[N] =
                    static_cast<to_type>(
                        fusion::at_c<1>(operands).template at<N, 0>()
                    );
            }
    };

    template <typename Unit, typename ValueType>
    struct matrix_to_temp_vector_assign<units::quantity<Unit, ValueType> >
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type::value_type to_type;
                fusion::at_c<0>(operands)[N] =
                    static_cast<to_type>(
                        fusion::at_c<1>(operands).template at<N, 0>().value()
                    );
            }
    };

    struct temp_vector_to_matrix_assign
    {
        template <std::size_t N, typename Operands>
        static void call (Operands const & operands)
            {
                typedef typename operand_n<Operands, 0>::type to_matrix;
                typedef typename result_of::value_at_c<to_matrix, N, 0>::type to_type;
                assign_from_value_type<to_type>::call(
                    fusion::at_c<0>(operands).template at<N, 0>(),
                    fusion::at_c<1>(operands)[N]
                );
            }
    };

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_DETAIL_ITERATION_HPP
