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
#include <boost/units_blas/exception.hpp>
#include <boost/units_blas/traits.hpp>
#include <boost/units_blas/detail/has_identity.hpp>
#include <boost/units_blas/detail/inverse_type.hpp>
#include <boost/units_blas/detail/lu.hpp>
#include <boost/units_blas/detail/one_value.hpp>
#include <boost/units_blas/detail/value_type.hpp>
#include <boost/units_blas/detail/zero_value.hpp>

#include <boost/throw_exception.hpp>
#include <boost/units/cmath.hpp>


namespace boost { namespace units_blas {

    namespace detail {

        // foldl
        template <std::size_t I, typename Fn, typename State>
        constexpr auto foldl_impl (Fn, State state, type_sequence<>)
        { return state; }

        template <std::size_t I,
                  typename Fn,
                  typename State,
                  typename Head,
                  typename ...Tail>
        constexpr auto foldl_impl (Fn f,
                                   State state,
                                   type_sequence<Head, Tail...>)
        {
            return foldl_impl<I + 1>(
                f,
                f.template call<I, Head>(state),
                type_sequence<Tail...>{}
            );
        }

        template <typename Fn, typename State, typename Head, typename ...Tail>
        constexpr auto foldl (Fn f,
                              State state,
                              type_sequence<Head, Tail...> seq)
        { return foldl_impl<0>(f, state, seq); }


        // column indices and types
        template <typename Matrix,
                  std::size_t X,
                  std::size_t Incr,
                  std::size_t I,
                  std::size_t N>
        struct indices_and_types_impl
        {
            template <typename Indices, typename Types>
            static constexpr auto call (std::pair<Indices, Types>)
            {
                using indices = decltype(push_back<X>(Indices{}));
                using type = typename std::tuple_element<
                    X,
                    typename Matrix::value_types
                >::type;
                using types = decltype(push_back<type>(Types{}));
                return indices_and_types_impl<
                    Matrix,
                    X + Incr,
                    Incr,
                    I + 1,
                    N
                >::call(std::pair<indices, types>{});
            }
        };

        template <typename Matrix,
                  std::size_t X,
                  std::size_t Incr,
                  std::size_t N>
        struct indices_and_types_impl<Matrix, X, Incr, N, N>
        {
            template <typename Seqs>
            static constexpr auto call (Seqs seqs)
            { return seqs; }
        };

        template <typename Matrix, std::size_t R>
        constexpr auto row_indices_and_types ()
        {
            return indices_and_types_impl<
                Matrix,
                R * Matrix::num_columns,
                1,
                0,
                Matrix::num_columns
            >::call(std::pair<index_sequence<>, type_sequence<>>{});
        }

        template <typename Matrix, std::size_t C>
        constexpr auto column_indices_and_types ()
        {
            return indices_and_types_impl<
                Matrix,
                C,
                Matrix::num_columns,
                0,
                Matrix::num_rows
            >::call(std::pair<index_sequence<>, type_sequence<>>{});
        }


        // indexed iteration
        template <std::size_t I,
                  typename F,
                  std::size_t Head,
                  std::size_t ...Tail>
        struct iterate_indexed_impl;

        template <std::size_t I,
                  typename F,
                  std::size_t Head,
                  std::size_t ...Tail>
        struct iterate_indexed_impl
        {
            static void call (F f, index_sequence<Head, Tail...>)
            {
                f.template call<I, Head>();
                iterate_indexed_impl<I + 1, F, Tail...>::call(
                    f,
                    index_sequence<Tail...>{}
                );
            }
        };

        template <std::size_t I, typename F, std::size_t Head>
        struct iterate_indexed_impl<I, F, Head>
        {
            static void call (F f, index_sequence<Head>)
            { f.template call<I, Head>(); }
        };

        template <typename F, std::size_t ...I>
        void iterate_indexed (F f, index_sequence<I...> seq)
        { iterate_indexed_impl<0, F, I...>::call(f, seq); }


        // row/column tuples
        template <typename Tuple, typename Matrix>
        struct tuple_assign
        {
            template <std::size_t I, std::size_t J>
            void call ()
            { std::get<I>(lhs_) = tuple_access::get<J>(rhs_); }

            Tuple & lhs_;
            Matrix const & rhs_;
        };

        template <std::size_t R, typename Matrix>
        constexpr auto row_tuple (Matrix m)
        {
            auto seqs = indices_and_types_impl<
                Matrix,
                R * Matrix::num_columns,
                1,
                0,
                Matrix::num_columns
            >::call(std::pair<index_sequence<>, type_sequence<>>{});
            auto retval = tuple_from_types(seqs.second);
            iterate_indexed(tuple_assign<decltype(retval), Matrix>{retval, m},
                            seqs.first);
            return retval;
        }

        template <std::size_t C, typename Matrix>
        constexpr auto column_tuple (Matrix m)
        {
            auto seqs = indices_and_types_impl<
                Matrix,
                C,
                Matrix::num_columns,
                0,
                Matrix::num_rows
            >::call(std::pair<index_sequence<>, type_sequence<>>{});
            auto retval = tuple_from_types(seqs.second);
            iterate_indexed(tuple_assign<decltype(retval), Matrix>{retval, m},
                            seqs.first);
            return retval;
        }


        // transpose indices and types
        template <typename Matrix, std::size_t I, std::size_t N>
        struct transpose_indices_and_types_impl
        {
            template <typename Indices, typename Types>
            static constexpr auto call (std::pair<Indices, Types>)
            {
                constexpr std::size_t row = I / Matrix::num_columns;
                constexpr std::size_t column = I % Matrix::num_columns;
                constexpr std::size_t transpose_i =
                    column * Matrix::num_rows + row;
                using indices = decltype(push_back<transpose_i>(Indices{}));

                using type = typename std::tuple_element<
                    transpose_i,
                    typename Matrix::value_types
                >::type;
                using types = decltype(push_back<type>(Types{}));

                return transpose_indices_and_types_impl<Matrix, I + 1, N>::call(
                    std::pair<indices, types>{}
                );
            }
        };

        template <typename Matrix, std::size_t N>
        struct transpose_indices_and_types_impl<Matrix, N, N>
        {
            template <typename Seqs>
            static constexpr auto call (Seqs seqs)
            { return seqs; }
        };


        // tuple dot product
        template <typename Tuple1, typename Tuple2>
        struct tuple_dot_impl
        {
            template <std::size_t I, typename T, typename State>
            auto call (State prev)
            { return prev + std::get<I>(t1_) * std::get<I>(t2_); }

            Tuple1 t1_;
            Tuple2 t2_;
        };

        template <typename Tuple1, typename Head, typename ...Tail>
        auto tuple_dot (Tuple1 t1, std::tuple<Head, Tail...> t2)
        {
            static_assert(
                std::tuple_size<Tuple1>::value ==
                std::tuple_size<std::tuple<Head, Tail...>>::value,
                "tuple_dot() must be given tuples of the same length"
            );
            auto state = std::get<0>(t1) * std::get<0>(t2);
            using function_object =
                tuple_dot_impl<Tuple1, std::tuple<Head, Tail...>>;
            function_object f{t1, t2};
            return foldl_impl<1>(f, state, type_sequence<Tail...>{});
        }


        // matrix slice
        template <std::size_t Head, std::size_t ...Tail>
        constexpr auto head (index_sequence<Head, Tail...>)
        { return Head; }

        template <std::size_t Head, std::size_t ...Tail>
        constexpr auto tail (index_sequence<Head, Tail...>)
        { return index_sequence<Tail...>{}; }

        template <typename Matrix,
                  std::size_t Row,
                  typename ColumnIndices,
                  std::size_t X>
        struct slice_indices_and_types_row_impl
        {
            template <typename Indices, typename Types>
            static constexpr auto call (std::pair<Indices, Types>)
            {
                using indices = decltype(push_back<X>(Indices{}));
                using type = typename std::tuple_element<
                    Row * Matrix::num_columns + head(ColumnIndices{}),
                    typename Matrix::value_types
                >::type;
                using types = decltype(push_back<type>(Types{}));
                return slice_indices_and_types_row_impl<
                    Matrix,
                    Row,
                    decltype(tail(ColumnIndices{})),
                    X + 1
                >::call(std::pair<indices, types>{});
            }
        };

        template <typename Matrix,
                  std::size_t Row,
                  std::size_t X>
        struct slice_indices_and_types_row_impl<
            Matrix,
            Row,
            index_sequence<>,
            X
        > {
            template <typename Seqs>
            static constexpr auto call (Seqs seqs)
            { return seqs; }
        };


        template <typename Matrix,
                  typename RowIndices,
                  typename ColumnIndices,
                  std::size_t X>
        struct slice_indices_and_types_impl
        {
            template <typename Indices, typename Types>
            static constexpr auto call (std::pair<Indices, Types> seqs)
            {
                auto new_seqs = slice_indices_and_types_row_impl<
                    Matrix,
                    head(RowIndices{}),
                    ColumnIndices,
                    X
                >::call(seqs);
                return slice_indices_and_types_impl<
                    Matrix,
                    decltype(tail(RowIndices{})),
                    ColumnIndices,
                    X + ColumnIndices::size
                >::call(new_seqs);
            }
        };

        template <typename Matrix,
                  typename ColumnIndices,
                  std::size_t X>
        struct slice_indices_and_types_impl<
            Matrix,
            index_sequence<>,
            ColumnIndices,
            X
        > {
            template <typename Seqs>
            static constexpr auto call (Seqs seqs)
            { return seqs; }
        };


        // matrix product types
        template <typename MatrixLHS,
                  typename MatrixRHS,
                  std::size_t I,
                  std::size_t N>
        struct matrix_product_types_impl
        {
            template <typename Types>
            static constexpr auto call (Types)
            {
                constexpr std::size_t row = I / MatrixRHS::num_columns;
                constexpr std::size_t column = I % MatrixRHS::num_columns;

                using type = decltype(
                    tuple_dot(
                        row_tuple<row>(std::declval<MatrixLHS>()),
                        column_tuple<column>(std::declval<MatrixRHS>())
                    )
                );
                using types = decltype(push_back<type>(Types{}));

                return matrix_product_types_impl<
                    MatrixLHS,
                    MatrixRHS,
                    I + 1,
                    N
                >::call(types{});
            }
        };

        template <typename MatrixLHS, typename MatrixRHS, std::size_t N>
        struct matrix_product_types_impl<MatrixLHS, MatrixRHS, N, N>
        {
            template <typename Seq>
            static constexpr auto call (Seq seq)
            { return seq; }
        };


        // elementwise product types
        template <typename MatrixLHS,
                  typename MatrixRHS,
                  std::size_t I,
                  std::size_t N>
        struct element_product_types_impl
        {
            template <typename Types>
            static constexpr auto call (Types)
            {
                using type = decltype(
                    tuple_access::get<I>(std::declval<MatrixLHS>()) *
                    tuple_access::get<I>(std::declval<MatrixRHS>())
                );
                using types = decltype(push_back<type>(Types{}));

                return element_product_types_impl<
                    MatrixLHS,
                    MatrixRHS,
                    I + 1,
                    N
                >::call(types{});
            }
        };

        template <typename MatrixLHS, typename MatrixRHS, std::size_t N>
        struct element_product_types_impl<MatrixLHS, MatrixRHS, N, N>
        {
            template <typename Seq>
            static constexpr auto call (Seq seq)
            { return seq; }
        };


        // elementwise product types
        template <typename MatrixLHS,
                  typename MatrixRHS,
                  std::size_t I,
                  std::size_t N>
        struct element_division_types_impl
        {
            template <typename Types>
            static constexpr auto call (Types)
            {
                using type = decltype(
                    tuple_access::get<I>(std::declval<MatrixLHS>()) /
                    tuple_access::get<I>(std::declval<MatrixRHS>())
                );
                using types = decltype(push_back<type>(Types{}));

                return element_division_types_impl<
                    MatrixLHS,
                    MatrixRHS,
                    I + 1,
                    N
                >::call(types{});
            }
        };

        template <typename MatrixLHS, typename MatrixRHS, std::size_t N>
        struct element_division_types_impl<MatrixLHS, MatrixRHS, N, N>
        {
            template <typename Seq>
            static constexpr auto call (Seq seq)
            { return seq; }
        };


        template <typename Matrix>
        struct swap
        {
            template <std::size_t I>
            void call ()
            {
                using std::swap;
                swap(tuple_access::get<I>(lhs_), tuple_access::get<I>(rhs_));
            }

            Matrix & lhs_;
            Matrix & rhs_;
        };

        template <typename Matrix>
        struct negate
        {
            template <std::size_t I>
            void call ()
            { tuple_access::get<I>(m_) = -tuple_access::get<I>(m_); }

            Matrix & m_;
        };

        template <typename ResultMatrix, typename MatrixLHS, typename MatrixRHS>
        struct prod
        {
            template <std::size_t I>
            void call ()
            {
                constexpr std::size_t row = I / ResultMatrix::num_columns;
                constexpr std::size_t column = I % ResultMatrix::num_columns;
                tuple_access::get<I>(m_) = tuple_dot(row_tuple<row>(lhs_),
                                                     column_tuple<column>(rhs_));
            }

            ResultMatrix & m_;
            MatrixLHS lhs_;
            MatrixRHS rhs_;
        };

        template <typename ResultMatrix, typename MatrixLHS, typename MatrixRHS>
        struct element_prod
        {
            template <std::size_t I>
            void call ()
            {
                tuple_access::get<I>(m_) =
                    tuple_access::get<I>(lhs_) * tuple_access::get<I>(rhs_);
            }

            ResultMatrix & m_;
            MatrixLHS lhs_;
            MatrixRHS rhs_;
        };

        template <typename ResultMatrix, typename MatrixLHS, typename MatrixRHS>
        struct element_div
        {
            template <std::size_t I>
            void call ()
            {
                tuple_access::get<I>(m_) =
                    tuple_access::get<I>(lhs_) / tuple_access::get<I>(rhs_);
            }

            ResultMatrix & m_;
            MatrixLHS lhs_;
            MatrixRHS rhs_;
        };

        template <typename Tuple, bool Abs>
        struct sum_impl
        {
            template <std::size_t I, typename T, typename State>
            auto call (State prev)
            {
                auto value = std::get<I>(t_);
                using std::abs;
                if (Abs)
                    value = abs(value);
                return prev + value;
            }

            Tuple t_;
        };

        template <typename Tuple>
        struct norm_2_impl
        {
            template <std::size_t I, typename T, typename State>
            auto call (State prev)
            {
                auto value = std::get<I>(t_);
                return prev + value * value;
            }

            Tuple t_;
        };

        template <typename Tuple>
        struct norm_inf_index_impl
        {
            template <std::size_t I, typename T, typename State>
            auto call (State prev)
            {
                using common_type = decltype(prev.second + std::get<I>(t_));
                using retval_type = std::pair<std::size_t, common_type>;
                using std::abs;
                auto value = abs(std::get<I>(t_));
                return
                    prev.second < value ?
                    retval_type{I, value} :
                    retval_type{prev.first, prev.second};
            }

            Tuple t_;
        };


        // matrix determinant type
#if 1 // TODO BOOST_UNITS_BLAS_USE_INEXACT_DETERMINANT_TYPE
        template <typename Matrix, std::size_t I, std::size_t N>
        struct simplified_determinant_type_impl
        {
            template <typename Prev>
            static constexpr auto call (Prev prev)
            {
                using type = typename std::tuple_element<
                    I * Matrix::num_columns + I,
                    typename Matrix::value_types
                >::type;
                return simplified_determinant_type_impl<Matrix, I + 1, N>::call(
                    prev * type{}
                );
            }
        };

        template <typename Matrix, std::size_t N>
        struct simplified_determinant_type_impl<Matrix, N, N>
        {
            template <typename Result>
            static constexpr auto call (Result result)
            { return result; }
        };

        template <typename Matrix>
        struct determinant_type
        {
            using first = typename std::tuple_element<
                0,
                typename Matrix::value_types
            >::type;
            using type = decltype(
                simplified_determinant_type_impl<
                    Matrix,
                    1,
                    Matrix::num_columns
                >::call(first{})
            );
        };
#else
        template <typename Matrix>
        struct determinant_type
        {
            using type = double; // TODO
        };
#endif


        // get value
        template <typename T, typename ValueType>
        struct get_value
        {
            static ValueType call (T t)
            { return static_cast<ValueType>(t); }
        };

        template <typename Unit, typename T, typename ValueType>
        struct get_value<units::quantity<Unit, T>, ValueType>
        {
            static ValueType call (units::quantity<Unit, T> u)
            { return static_cast<ValueType>(u.value()); }
        };


        // matrix_t <--> temp matrix/array
        template <typename TempMatrix, typename Matrix>
        struct assign_to_temp_matrix
        {
            template <std::size_t I>
            void call ()
            {
                constexpr std::size_t row = I / Matrix::num_columns;
                constexpr std::size_t column = I % Matrix::num_columns;
                tmp_[row][column] = get_value<
                    typename std::tuple_element<
                        I,
                        typename Matrix::value_types
                    >::type,
                    typename TempMatrix::value_type::value_type
                >::call(tuple_access::get<I>(m_));
            }

            TempMatrix & tmp_;
            Matrix m_;
        };

        template <typename T, typename ValueType>
        struct make_value
        {
            static auto call (ValueType value)
            { return static_cast<T>(value); }
        };

        template <typename Unit, typename T, typename ValueType>
        struct make_value<units::quantity<Unit, T>, ValueType>
        {
            static auto call (ValueType value)
            {
                return units::quantity<Unit, T>::from_value(
                    static_cast<T>(value)
                );
            }
        };

        template <typename TempMatrix, typename Matrix>
        struct assign_from_temp_matrix
        {
            template <std::size_t I>
            void call ()
            {
                constexpr std::size_t row = I / Matrix::num_columns;
                constexpr std::size_t column = I % Matrix::num_columns;
                tuple_access::get<I>(m_) = make_value<
                    typename std::tuple_element<
                        I,
                        typename Matrix::value_types
                    >::type,
                    typename TempMatrix::value_type::value_type
                >::call(tmp_[row][column]);
            }

            TempMatrix tmp_;
            Matrix & m_;
        };


        template <typename TempArray, typename Matrix>
        struct assign_to_temp_array
        {
            template <std::size_t I>
            void call ()
            {
                tmp_[I] = get_value<
                    typename std::tuple_element<
                        I,
                        typename Matrix::value_types
                    >::type,
                    typename TempArray::value_type
                >::call(tuple_access::get<I>(m_));
            }

            TempArray & tmp_;
            Matrix m_;
        };

        template <typename TempArray, typename Matrix>
        struct assign_from_temp_array
        {
            template <std::size_t I>
            void call ()
            {
                tuple_access::get<I>(m_) = make_value<
                    typename std::tuple_element<
                        I,
                        typename Matrix::value_types
                    >::type,
                    typename TempArray::value_type
                >::call(tmp_[I]);
            }

            TempArray tmp_;
            Matrix & m_;
        };


        // inverse() support
        template <typename TempMatrix, typename Indices, std::size_t Columns>
        struct assign_inverted_column
        {
            template <std::size_t C>
            void call ()
            {
                using value_type = typename TempMatrix::value_type::value_type;
                std::array<value_type, Columns> columns = {};
                columns[C] = value_type{1.0};
                detail::lu_substitute(from_, indices_, columns);
                for (std::size_t i = 0; i < Columns; ++i) {
                    to_[i][C] = columns[i];
                }
            }

            TempMatrix & to_;
            TempMatrix from_;
            Indices indices_;
        };

    }

    /** Returns a @c matrix<> consisting of only the rows and columns of @c m
        specified by @c RowIndices and @c ColumnIndices.  @c m must be a @c
        matrix<>.  Both of @c RowIndices and @c ColumnIndices must be an @c
        index_sequence<> ; all values in @c RowIndices and @c ColumnIndices must
        be less than the number of rows and columns in @c m, respectively.  Note
        that duplication and order preservation are not enforced for the values
        in @c RowIndices and @c ColumnIndices.  It is therefore possible to use
        @c slice<>() to rearrange and/or duplicate rows and/or columns. */
    template <typename RowIndices,
              typename ColumnIndices,
              typename Tuple,
              std::size_t Rows,
              std::size_t Columns>
    auto slice (matrix_t<Tuple, Rows, Columns> m)
    {
        using matrix_type = matrix_t<Tuple, Rows, Columns>;

        static_assert(
            0 < RowIndices::size,
            "slice() requires at least one row index to compute its result"
        );

        static_assert(
            0 < ColumnIndices::size,
            "slice() requires at least one column index to compute its result"
        );

        auto seqs = detail::slice_indices_and_types_impl<
            matrix_type,
            RowIndices,
            ColumnIndices,
            0
        >::call(std::pair<index_sequence<>, detail::type_sequence<>>{});

        auto tuple = detail::tuple_from_types(seqs.second);
        iterate_indexed(
            detail::tuple_assign<decltype(tuple), matrix_type>{tuple, m},
            seqs.first
        );

        return
            detail::make_matrix<RowIndices::size, ColumnIndices::size>(tuple);
    }

    /** Returns a const-preserved reference to the element at row @c I, column
        @c J of @c m.  @c m must be a @c matrix<>. */
    template <std::size_t I, std::size_t J, typename Matrix>
    decltype(auto) at (Matrix && m)
    { return m.template at<I, J>(); }

    /** Returns the tranpose of @c m.  @c m must be a @c matrix<>. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    auto transpose (matrix_t<Tuple, Rows, Columns> m)
    {
        using matrix_type = matrix_t<Tuple, Rows, Columns>;
        auto seqs = detail::transpose_indices_and_types_impl<
            matrix_type,
            0,
            matrix_type::num_elements
        >::call(std::pair<index_sequence<>, detail::type_sequence<>>{});
        auto tuple = tuple_from_types(seqs.second);
        iterate_indexed(
            detail::tuple_assign<decltype(tuple), matrix_type>{tuple, m},
            seqs.first
        );
        return detail::make_matrix<Columns, Rows>(tuple);
    }

    /** Returns the negation of @c m.  @c m must be a @c matrix<>. */
    template <std::size_t Rows, std::size_t Columns, typename ...T>
    auto neg (matrix_t<std::tuple<T...>, Rows, Columns> m)
    {
        detail::iterate_simple<m.num_elements>(
            detail::negate<matrix_t<std::tuple<T...>, Rows, Columns>>{m}
        );
        return m;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the negation of @c m.  @c m must be a @c matrix<>. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    auto operator- (matrix_t<Tuple, Rows, Columns> m) -> decltype(neg(m))
    { return neg(m); }

#endif

    /** Returns the elementwise sum of @c lhs and @c rhs. @c lhs and @c rhs
        must be <c>matrix<></c>s with the same dimensions.  Also, every sum
        <c>lhs(i, j) + rhs(i, j)</c> must be a valid operation.  */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns>
    auto sum (matrix_t<Tuple1, Rows, Columns> lhs,
              matrix_t<Tuple2, Rows, Columns> rhs)
    { return lhs += rhs; }

    /** Returns the elementwise difference of @c lhs and @c rhs. @c lhs and @c
        rhs must be <c>matrix<></c>s with the same dimensions.  Also, every
        difference <c>lhs(i, j) - rhs(i, j)</c> must be a valid operation.  */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns>
    auto diff (matrix_t<Tuple1, Rows, Columns> lhs,
               matrix_t<Tuple2, Rows, Columns> rhs)
    { return lhs -= rhs; }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the elementwise sum of @c lhs and @c rhs. @c lhs and @c rhs
        must be <c>matrix<></c>s with the same dimensions.  Also, every sum
        <c>lhs(i, j) + rhs(i, j)</c> must be a valid operation.  */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns>
    auto operator+ (matrix_t<Tuple1, Rows, Columns> lhs,
                    matrix_t<Tuple2, Rows, Columns> rhs) -> decltype(sum(lhs, rhs))
    { return sum(lhs, rhs); }

    /** Returns the elementwise difference of @c lhs and @c rhs. @c lhs and @c
        rhs must be <c>matrix<></c>s with the same dimensions.  Also, every
        difference <c>lhs(i, j) - rhs(i, j)</c> must be a valid operation.  */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns>
    auto operator- (matrix_t<Tuple1, Rows, Columns> lhs,
                    matrix_t<Tuple2, Rows, Columns> rhs) -> decltype(diff(lhs, rhs))
    { return diff(lhs, rhs); }

#endif

    /** Returns the matrix-product of @c lhs and @c rhs. @c lhs and @c rhs
        must be <c>matrix<></c>s, and the number of columns in @c lhs must the
        same as the number of rows in @c rhs.  Also, a matrix-product type
        must exist for @c lhs and @c rhs (some otherwise-suitable pairs of
        <c>matrix<></c>s do not have a matrix-product that makes sense when
        their elements are unit types). */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows1,
              std::size_t Inner,
              std::size_t Columns2>
    auto prod (matrix_t<Tuple1, Rows1, Inner> lhs,
               matrix_t<Tuple2, Inner, Columns2> rhs)
    {
        using lhs_type = matrix_t<Tuple1, Rows1, Inner>;
        using rhs_type = matrix_t<Tuple2, Inner, Columns2>;
        auto seq = detail::matrix_product_types_impl<
            lhs_type,
            rhs_type,
            0,
            Rows1 * Columns2
        >::call(detail::type_sequence<>{});
        auto retval =
            detail::make_matrix<Rows1, Columns2>(tuple_from_types(seq));
        detail::iterate_simple<Rows1 * Columns2>(
            detail::prod<decltype(retval), lhs_type, rhs_type>{retval, lhs, rhs}
        );
        return retval;
    }

    /** Returns the product of @c m and @c t.  @c m must be a @c matrix<>, and
        @c T must not be a @c matrix<>. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns, typename T>
    auto prod (matrix_t<Tuple, Rows, Columns> m,
               T t,
               typename std::enable_if<
                   !is_matrix_or_derived<T>::value
               >::type* = 0)
    {
        matrix_t<Tuple, Rows, Columns> retval = m;
        retval *= t;
        return retval;
    }

    /** Returns the product of @c t and @c m.  @c m must be a @c matrix<>, and
        @c T must not be a @c matrix<>. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns, typename T>
    auto prod (T t,
               matrix_t<Tuple, Rows, Columns> m,
               typename std::enable_if<
                   !is_matrix_or_derived<T>::value
               >::type* = 0)
    { return prod(m, t); }

    /** Returns the result of dividing @c m by @c t.  @c m must be a @c
        matrix<>, and @c T must not be a @c matrix<>. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns, typename T>
    auto div (matrix_t<Tuple, Rows, Columns> m, T t)
    {
        matrix_t<Tuple, Rows, Columns> retval = m;
        retval /= t;
        return retval;
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the matrix-product of @c lhs and @c rhs. @c lhs and @c rhs
        must be <c>matrix<></c>s, and the number of columns in @c lhs must the
        same as the number of rows in @c rhs.  Also, a matrix-product type
        must exist for @c lhs and @c rhs (some otherwise-suitable pairs of
        <c>matrix<></c>s do not have a matrix-product that makes sense when
        their elements are unit types). */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows1,
              std::size_t Inner,
              std::size_t Columns2>
    auto operator* (matrix_t<Tuple1, Rows1, Inner> lhs,
                    matrix_t<Tuple2, Inner, Columns2> rhs) -> decltype(prod(lhs, rhs))
    { return prod(lhs, rhs); }

    /** Returns the product of @c m and @c t.  @c m must be a @c matrix<>, and
        @c T must not be a @c matrix<>. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns, typename T,
              typename Enable = typename std::enable_if<
                  !is_matrix_or_derived<T>::value
              >::type>
    auto operator* (matrix_t<Tuple, Rows, Columns> m, T t) -> decltype(prod(m, t))
    { return prod(m, t); }

    /** Returns the product of @c t and @c m.  @c m must be a @c matrix<>, and
        @c T must not be a @c matrix<>. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns, typename T,
              typename Enable = typename std::enable_if<
                  !is_matrix_or_derived<T>::value
              >::type>
    auto operator* (T t, matrix_t<Tuple, Rows, Columns> m) -> decltype(prod(m, t))
    { return prod(m, t); }

    /** Returns the result of dividing @c m by @c t.  @c m must be a @c
        matrix<>, and @c T must not be a @c matrix<>. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns, typename T>
    auto operator/ (matrix_t<Tuple, Rows, Columns> m, T t) -> decltype(div(m, t))
    { return div(m, t); }

#endif

    /** Returns the elementwise multiplication of the elements of @c lhs by
        the elements of @c rhs.  @c lhs and @c rhs must be <c>matrix<></c>s
        with the same dimensions. */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns>
    auto element_prod (matrix_t<Tuple1, Rows, Columns> lhs,
                       matrix_t<Tuple2, Rows, Columns> rhs)
    {
        using lhs_type = matrix_t<Tuple1, Rows, Columns>;
        using rhs_type = matrix_t<Tuple2, Rows, Columns>;
        auto seq = detail::element_product_types_impl<
            lhs_type,
            rhs_type,
            0,
            Rows * Columns
        >::call(detail::type_sequence<>{});
        auto retval =
            detail::make_matrix<Rows, Columns>(tuple_from_types(seq));
        detail::iterate_simple<Rows * Columns>(
            detail::element_prod<
                decltype(retval),
                lhs_type,
                rhs_type
            >{retval, lhs, rhs}
        );
        return retval;
    }

    /** Returns the elementwise division of the elements of @c lhs by the
        elements of @c rhs.  @c lhs and @c rhs must be <c>matrix<></c>s with
        the same dimensions. */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns>
    auto element_div (matrix_t<Tuple1, Rows, Columns> lhs,
                      matrix_t<Tuple2, Rows, Columns> rhs)
    {
        using lhs_type = matrix_t<Tuple1, Rows, Columns>;
        using rhs_type = matrix_t<Tuple2, Rows, Columns>;
        auto seq = detail::element_division_types_impl<
            lhs_type,
            rhs_type,
            0,
            Rows * Columns
        >::call(detail::type_sequence<>{});
        auto retval =
            detail::make_matrix<Rows, Columns>(tuple_from_types(seq));
        detail::iterate_simple<Rows * Columns>(
            detail::element_div<
                decltype(retval),
                lhs_type,
                rhs_type
            >{retval, lhs, rhs}
        );
        return retval;
    }

    /** Swaps the values in @c lhs and @c rhs.  Note that this is an
        O(@c size<matrix<T> >::value) operation. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    void swap (matrix_t<Tuple, Rows, Columns> & lhs,
               matrix_t<Tuple, Rows, Columns> & rhs)
    {
        using matrix_type = matrix_t<Tuple, Rows, Columns>;
        detail::iterate_simple<matrix_type::num_elements>(
            detail::swap<matrix_type>{lhs, rhs}
        );
    }

    /** Returns the dot product of @c lhs and @c rhs.  @c VectorL and @c VectorR
        must be <c>matrix<></c>s, and must have the same dimensions.
        Additionally, both <c>matrix<></c>s must be "vectors" of length greater
        than 1. */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns>
    auto dot (matrix_t<Tuple1, Rows, Columns> lhs,
              matrix_t<Tuple2, Rows, Columns> rhs,
              typename std::enable_if<
                  Columns == 1 && 1 < Rows || Rows == 1 && 1 < Columns
              >::type* = 0)
    {
        return detail::tuple_dot(detail::tuple_access::all(lhs),
                                 detail::tuple_access::all(rhs));
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the dot product of @c lhs and @c rhs.  @c VectorL and @c VectorR
        must be <c>matrix<></c>s, and must have the same dimensions.
        Additionally, both <c>matrix<></c>s must be "vectors" of length greater
        than 1. */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns,
              typename Enable = typename std::enable_if<
                  Columns == 1 && 1 < Rows || Rows == 1 && 1 < Columns
              >::type>
    auto operator* (matrix_t<Tuple1, Rows, Columns> lhs,
                    matrix_t<Tuple2, Rows, Columns> rhs) -> decltype(dot(lhs, rhs))
    { return dot(lhs, rhs); }

#endif

    /** Returns the cross product of @c lhs with @c rhs.  @c VectorL and @c
        VectorR must both be <c>matrix<></c>s, and must both be 3 x 1
        "vectors".  Also, a cross product type must exist for @c VectorL and
        @c VectorR (some otherwise-suitable pairs of <c>matrix<></c>s do not
        have a cross product that makes sense when their elements are unit
        types). */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns>
    auto cross (matrix_t<Tuple1, Rows, Columns> lhs,
                matrix_t<Tuple2, Rows, Columns> rhs,
                typename std::enable_if<
                    Columns == 3 && Rows == 1 || Rows == 3 && Columns == 1
                >::type* = 0)
    {
        auto l = detail::tuple_access::all(lhs);
        auto r = detail::tuple_access::all(rhs);
        auto _0 =
            std::get<1>(l) * std::get<2>(r) -
            std::get<2>(l) * std::get<1>(r);
        auto _1 =
            std::get<2>(l) * std::get<0>(r) -
            std::get<0>(l) * std::get<2>(r);
        auto _2 =
            std::get<0>(l) * std::get<1>(r) -
            std::get<1>(l) * std::get<0>(r);
        return detail::make_matrix<Rows, Columns>(std::make_tuple(_0, _1, _2));
    }

#if BOOST_UNITS_BLAS_USE_OPERATORS_FOR_MATRIX_OPERATIONS

    /** Returns the cross product of @c lhs with @c rhs.  @c VectorL and @c
        VectorR must both be <c>matrix<></c>s, and must both be 3 x 1
        "vectors".  Also, a cross product type must exist for @c VectorL and
        @c VectorR (some otherwise-suitable pairs of <c>matrix<></c>s do not
        have a cross product that makes sense when their elements are unit
        types). */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns,
              typename Enable = typename std::enable_if<
                  Columns == 3 && Rows == 1 || Rows == 3 && Columns == 1
              >::type>
    auto operator^ (matrix_t<Tuple1, Rows, Columns> lhs,
                    matrix_t<Tuple2, Rows, Columns> rhs) -> decltype(cross(lhs, rhs))
    { return cross(lhs, rhs); }

#endif

    /** Returns the sum of all elements in @c v.  @c Vector must be a "vector"
        @c matrix<>.  Also, a sum type must exist for @c Vector (some
        otherwise-suitable <c>matrix<></c>s do not have a sum that makes sense
        when their elements are unit types). */
    template <std::size_t Rows,
              std::size_t Columns,
              typename Head,
              typename ...Tail>
    auto sum (matrix_t<std::tuple<Head, Tail...>, Rows, Columns> v,
              typename std::enable_if<
                  Rows == 1 || Columns == 1
              >::type* = 0)
    {
        auto f = detail::sum_impl<std::tuple<Head, Tail...>, false>{
            detail::tuple_access::all(v)
        };
        auto state = detail::tuple_access::get<0>(v);
        return detail::foldl_impl<1>(
            f,
            state,
            detail::type_sequence<Tail...>{}
        );
    }

    /** Returns the sum of the absolute values of all elements in @c v.  @c
        Vector must be a "vector" @c matrix<>.  Also, a sum type must exist
        for @c Vector (some otherwise-suitable <c>matrix<></c>s do not have a
        sum that makes sense when their elements are unit types). */
    template <std::size_t Rows,
              std::size_t Columns,
              typename Head,
              typename ...Tail>
    auto norm_1 (matrix_t<std::tuple<Head, Tail...>, Rows, Columns> v,
                 typename std::enable_if<
                     Rows == 1 || Columns == 1
                 >::type* = 0)
    {
        auto f = detail::sum_impl<std::tuple<Head, Tail...>, true>{
            detail::tuple_access::all(v)
        };
        using std::abs;
        auto state = abs(detail::tuple_access::get<0>(v));
        return detail::foldl_impl<1>(
            f,
            state,
            detail::type_sequence<Tail...>{}
        );
    }

    /** Returns the square root of the sum of the squares of all elements in
        @c v.  @c Vector must be a "vector" @c matrix<>.  Also, a sum type
        must exist for @c Vector (some otherwise-suitable <c>matrix<></c>s do
        not have a sum that makes sense when their elements are unit
        types). */
    template <std::size_t Rows,
              std::size_t Columns,
              typename Head,
              typename ...Tail>
    auto norm_2 (matrix_t<std::tuple<Head, Tail...>, Rows, Columns> v,
                 typename std::enable_if<
                     Rows == 1 || Columns == 1
                 >::type* = 0)
    {
        auto f = detail::norm_2_impl<std::tuple<Head, Tail...>>{
            detail::tuple_access::all(v)
        };
        auto first = detail::tuple_access::get<0>(v);
        auto state = first * first;
        auto sum = detail::foldl_impl<1>(
            f,
            state,
            detail::type_sequence<Tail...>{}
        );
        using std::sqrt;
        return sqrt(sum);
    }

    /** Returns the max of the absolute values of all elements in @c v.  @c
        Vector must be a "vector" @c matrix<>.  Also, @c operator+ and @c
        operator< must be defined for all pairs of elements in @c Vector. */
    template <std::size_t Rows, std::size_t Columns, typename ...T>
    auto norm_inf (matrix_t<std::tuple<T...>, Rows, Columns> v,
                   typename std::enable_if<
                       Rows == 1 || Columns == 1
                   >::type* = 0)
    {
        auto f = detail::norm_inf_index_impl<std::tuple<T...>>{
            detail::tuple_access::all(v)
        };
        using std::abs;
        auto state = std::make_pair(
            std::size_t{0},
            abs(detail::tuple_access::get<0>(v))
        );
        return detail::foldl(f, state, detail::type_sequence<T...>{}).second;
    }

    /** Returns the index of the first element in @c v equal to @c
        norm_inf(v).  @c Vector must be a "vector" @c matrix<>.  Also, @c
        operator+ and @c operator< must be defined for all pairs of elements
        in @c Vector. */
    template <std::size_t Rows, std::size_t Columns, typename ...T>
    auto norm_inf_index (matrix_t<std::tuple<T...>, Rows, Columns> v,
                         typename std::enable_if<
                             Rows == 1 || Columns == 1
                         >::type* = 0)
    {
        auto f = detail::norm_inf_index_impl<std::tuple<T...>>{
            detail::tuple_access::all(v)
        };
        using std::abs;
        auto state = std::make_pair(
            std::size_t{0},
            abs(detail::tuple_access::get<0>(v))
        );
        return detail::foldl(f, state, detail::type_sequence<T...>{}).first;
    }

#ifndef BOOST_UNITS_BLAS_DOXYGEN
    template <typename Tuple>
    auto determinant (matrix_t<Tuple, 1, 1> m)
    { return m.template at<0, 0>(); }

    template <typename Tuple>
    auto determinant (matrix_t<Tuple, 2, 2> m)
    {
        return
            m.template at<0, 0>() * m.template at<1, 1>() -
            m.template at<0, 1>() * m.template at<1, 0>();
    }

    template <typename Tuple>
    auto determinant (matrix_t<Tuple, 3, 3> m)
    {
        return
            m.template at<0, 0>() *
            (m.template at<1, 1>() * m.template at<2, 2>() -
             m.template at<1, 2>() * m.template at<2, 1>()) -
            m.template at<0, 1>() *
            (m.template at<1, 0>() * m.template at<2, 2>() -
             m.template at<1, 2>() * m.template at<2, 0>()) +
            m.template at<0, 2>() *
            (m.template at<1, 0>() * m.template at<2, 1>() -
             m.template at<1, 1>() * m.template at<2, 0>());
    }
#endif

    /** Returns the determinant of @c m.  @c m must be a @c matrix<>, and must
        be square.  Also, a determinant type must exist for @c m (some
        otherwise-suitable <c>matrix<></c>s do not have a determinant that
        makes sense when their elements are unit types).  */
    template <typename Tuple, std::size_t Rows>
    auto determinant (matrix_t<Tuple, Rows, Rows> m,
                      typename std::enable_if<
                          4 <= Rows
                      >::type* = 0)
    {
        using matrix_type = matrix_t<Tuple, Rows, Rows>;

        static_assert(
            detail::has_uniform_dimensional_units<matrix_type>::value,
            "LU decomposition requires a matrix that does not have mixed "
            "units of the same dimension (e.g. centimeters and meters in the "
            "same matrix won't work)."
        );

        using result_type =
            typename detail::determinant_type<matrix_type>::type;
        using raw_value_type = typename detail::value_type<result_type>::type;
        using temp_matrix_type =
            std::array<std::array<raw_value_type, Rows>, Rows>;

        temp_matrix_type temp_matrix;
        detail::iterate_simple<Rows * Rows>(
            detail::assign_to_temp_matrix<
                temp_matrix_type,
                matrix_type
            >{temp_matrix, m}
        );

        result_type retval = detail::zero_value<result_type>::value();

        std::array<std::size_t, Rows> indices;
        auto lu_result = detail::lu_decompose(temp_matrix, indices);
        if (lu_result.second) {
            raw_value_type tmp = lu_result.first;
            for (std::size_t i = 0; i < Rows; ++i) {
                tmp *= temp_matrix[i][i];
            }
            retval = detail::one_value<result_type>::value() * tmp;
        }

        return retval;
    }

    /** Returns the inverse of @c m.  Throws @c singular_matrix if @c m is
        found to be singular.  @c m must be a @c matrix<>, and must be
        square. */
    template <typename Tuple, std::size_t Rows>
    auto inverse (matrix_t<Tuple, Rows, Rows> m)
    {
        using matrix_type = matrix_t<Tuple, Rows, Rows>;

        static_assert(
            detail::has_identity<matrix_type>::value,
            "The given matrix type has no valid inverse, because it has no "
            "identity-matrix type.   A valid inverse of matrix type M must "
            "yield an identity matrix I such that M = M * I is a valid "
            "operation."
        );

        static_assert(
            detail::has_uniform_dimensional_units<matrix_type>::value,
            "LU decomposition requires a matrix that does not have mixed "
            "units of the same dimension (e.g. centimeters and meters in the "
            "same matrix won't work)."
        );

        using result_type = typename detail::inverse_type<matrix_type>::type;
        using temp_value_type =
            typename detail::determinant_type<matrix_type>::type;
        using raw_value_type =
            typename detail::value_type<temp_value_type>::type;
        using temp_matrix_type =
            std::array<std::array<raw_value_type, Rows>, Rows>;

        temp_matrix_type temp_matrix;
        detail::iterate_simple<Rows * Rows>(
            detail::assign_to_temp_matrix<
                temp_matrix_type,
                matrix_type
            >{temp_matrix, m}
        );

        result_type retval;

        std::array<std::size_t, Rows> indices;
        bool nonsingular = detail::lu_decompose(temp_matrix, indices).second;
        if (!nonsingular)
            throw_exception(singular_matrix{});

        temp_matrix_type temp_result_matrix;
        detail::iterate_simple<Rows>(
            detail::assign_inverted_column<
                temp_matrix_type,
                std::array<std::size_t, Rows>,
                Rows
            >{temp_result_matrix, temp_matrix, indices}
        );
        detail::iterate_simple<Rows * Rows>(
            detail::assign_from_temp_matrix<
                temp_matrix_type,
                result_type
            >{temp_result_matrix, retval}
        );

        return retval;
    }

    /** Returns the solution to the equation Ax = b in @c x.  Throws @c
        singular_matrix if @c A is found to be singular.  @c A must be a @c
        matrix<>, and must be square.  @c x and @c b must be "vector"
        <c>matrix<></c>s with the same dimensions, and must have a number of
        rows equal to the number of columns in @c A. */
    template <typename Tuple1,
              typename Tuple2,
              typename Tuple3,
              std::size_t Rows>
    void solve (matrix_t<Tuple1, Rows, Rows> A,
                matrix_t<Tuple2, Rows, 1> b,
                matrix_t<Tuple3, Rows, 1> & x)
    {
        using a_type = matrix_t<Tuple1, Rows, Rows>;
        using b_type = matrix_t<Tuple2, Rows, 1>;
        using x_type = matrix_t<Tuple3, Rows, 1>;

        static_assert(
            detail::has_uniform_dimensional_units<a_type>::value,
            "LU decomposition requires a matrix that does not have mixed "
            "units of the same dimension (e.g. centimeters and meters in the "
            "same matrix won't work)."
        );

        using temp_value_type = typename detail::determinant_type<a_type>::type;
        using raw_value_type = typename detail::value_type<temp_value_type>::type;

        using temp_matrix_type =
            std::array<std::array<raw_value_type, Rows>, Rows>;
        temp_matrix_type temp_A;
        detail::iterate_simple<Rows * Rows>(
            detail::assign_to_temp_matrix<
                temp_matrix_type,
                a_type
            >{temp_A, A}
        );

        std::array<std::size_t, Rows> indices;
        bool nonsingular = detail::lu_decompose(temp_A, indices).second;
        if (!nonsingular)
            throw_exception(singular_matrix{});

        using temp_vector_type = std::array<raw_value_type, Rows>;
        temp_vector_type temp_vector;
        detail::iterate_simple<Rows>(
            detail::assign_to_temp_array<
                temp_vector_type,
                b_type
            >{temp_vector, b}
        );

        detail::lu_substitute(temp_A, indices, temp_vector);

        detail::iterate_simple<Rows>(
            detail::assign_from_temp_array<
                temp_vector_type,
                x_type
            >{temp_vector, x}
        );
    }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_OPERATIONS_HPP
