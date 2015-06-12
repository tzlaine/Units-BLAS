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

    template <typename Tuple>
    auto determinant (matrix_t<Tuple, 3, 3> m);

    namespace detail {

        // compile-time foldl
        template <typename F, typename State, typename Seq>
        struct static_foldl;

        template <typename F, typename State, typename Head, typename ...Tail>
        struct static_foldl<F, State, type_sequence<Head, Tail...>>
        {
            using type = typename static_foldl<
                F,
                typename F::template apply<State, Head>::type,
                type_sequence<Tail...>
            >::type;
        };

        template <typename F, typename State>
        struct static_foldl<F, State, type_sequence<>>
        {
            using type = State;
        };

        template <typename F, typename State, typename Seq>
        using static_foldl_t = typename static_foldl<F, State, Seq>::type;


        // run-time foldl
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


        // column/row indices and types
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
                using type = tuple_element_t<X, Matrix>;
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
            static void call (F f, std::index_sequence<Head, Tail...>)
            {
                f.template call<I, Head>();
                iterate_indexed_impl<I + 1, F, Tail...>::call(
                    f,
                    std::index_sequence<Tail...>{}
                );
            }
        };

        template <std::size_t I, typename F, std::size_t Head>
        struct iterate_indexed_impl<I, F, Head>
        {
            static void call (F f, std::index_sequence<Head>)
            { f.template call<I, Head>(); }
        };

        template <typename F, std::size_t ...I>
        void iterate_indexed (F f, std::index_sequence<I...> seq)
        { iterate_indexed_impl<0, F, I...>::call(f, seq); }


        // row/column tuples
        template <typename Tuple, typename Matrix>
        struct tuple_assign
        {
            template <std::size_t I, std::size_t J>
            void call ()
            { hana::at_c<I>(lhs_) = tuple_access::get<J>(rhs_); }

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
            >::call(std::pair<std::index_sequence<>, type_sequence<>>{});
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
            >::call(std::pair<std::index_sequence<>, type_sequence<>>{});
            auto retval = tuple_from_types(seqs.second);
            iterate_indexed(tuple_assign<decltype(retval), Matrix>{retval, m},
                            seqs.first);
            return retval;
        }


        // transpose indices and types
        template <typename Matrix, std::size_t ...I>
        auto transpose_indices_and_types (std::index_sequence<I...>)
        {
            return std::pair<
                std::index_sequence<transpose_index<Matrix>(I)...>,
                type_sequence<
                    tuple_element_t<transpose_index<Matrix>(I), Matrix>...
                >
            >{};
        }


        // tuple dot product
        template <typename Tuple1, typename Tuple2>
        struct tuple_dot_impl
        {
            template <std::size_t I, typename T, typename State>
            auto call (State prev)
            { return prev + hana::at_c<I>(t1_) * hana::at_c<I>(t2_); }

            Tuple1 t1_;
            Tuple2 t2_;
        };

        template <typename Tuple1, typename Head, typename ...Tail>
        auto tuple_dot (Tuple1 t1, hana::_tuple<Head, Tail...> t2)
        {
            static_assert(
                Tuple1::size ==
                hana::_tuple<Head, Tail...>::size,
                "tuple_dot() must be given tuples of the same length"
            );
            auto state = hana::at_c<0>(t1) * hana::at_c<0>(t2);
            using function_object =
                tuple_dot_impl<Tuple1, hana::_tuple<Head, Tail...>>;
            function_object f{t1, t2};
            return foldl_impl<1>(f, state, type_sequence<Tail...>{});
        }


        // matrix slice
        template <typename Matrix,
                  std::size_t Row,
                  std::size_t ...ColumnIndices>
        struct slice_indices_and_types_row_impl
        {
            template <std::size_t ...I>
            static constexpr auto indices (std::index_sequence<I...>)
            {
                return std::index_sequence<
                    I...,
                    (Row * Matrix::num_columns + ColumnIndices)...
                >{};
            }

            template <typename ...T>
            static constexpr auto types (type_sequence<T...>)
            {
                return type_sequence<
                    T...,
                    tuple_element_t<
                        Row * Matrix::num_columns + ColumnIndices,
                        Matrix
                    >...
                >{};
            }
        };

        template <typename Matrix,
                  std::size_t Row,
                  typename Seqs,
                  std::size_t ...ColumnIndices>
        constexpr auto slice_indices_and_types_row (
            Seqs seqs,
            std::index_sequence<ColumnIndices...>
        ) {
            using impl = slice_indices_and_types_row_impl<
                Matrix,
                Row,
                ColumnIndices...
            >;
            return std::make_pair(
                impl::indices(seqs.first),
                impl::types(seqs.second)
            );
        }

        template <typename Matrix,
                  typename ColumnIndices,
                  typename Seqs,
                  std::size_t HeadRowIndex,
                  std::size_t ...TailRowIndices>
        constexpr auto slice_indices_and_types (
            Seqs seqs,
            std::index_sequence<HeadRowIndex, TailRowIndices...>
        ) {
            return slice_indices_and_types<Matrix, ColumnIndices>(
                slice_indices_and_types_row<Matrix, HeadRowIndex>(
                    seqs,
                    ColumnIndices{}
                ),
                std::index_sequence<TailRowIndices...>{}
            );
        }

        template <typename Matrix,
                  typename ColumnIndices,
                  typename Seqs>
        constexpr auto slice_indices_and_types (Seqs seqs,
                                                std::index_sequence<>)
        { return seqs; }


        // matrix product
        template <typename MatrixLHS, typename MatrixRHS, std::size_t ...I>
        auto matrix_prod_impl (MatrixLHS lhs,
                               MatrixRHS rhs,
                               std::index_sequence<I...>)
        {
            return hana::make<hana::Tuple>(
                tuple_dot(
                    row_tuple<I / MatrixRHS::num_columns>(lhs),
                    column_tuple<I % MatrixRHS::num_columns>(rhs)
                )...
            );
        }


        // scalar operations
        template <typename Matrix, typename T, std::size_t ...I>
        auto scalar_prod_impl (Matrix m, T t, std::index_sequence<I...>)
        { return hana::make<hana::Tuple>((tuple_access::get<I>(m) * t)...); }

        template <typename Matrix, typename T, std::size_t ...I>
        auto scalar_div_impl (Matrix m, T t, std::index_sequence<I...>)
        { return hana::make<hana::Tuple>((tuple_access::get<I>(m) / t)...); }


        // elementwise operations
        template <typename MatrixLHS, typename MatrixRHS, std::size_t ...I>
        auto element_sum_impl (MatrixLHS lhs,
                               MatrixRHS rhs,
                               std::index_sequence<I...>)
        {
            return hana::make<hana::Tuple>(
                (tuple_access::get<I>(lhs) + tuple_access::get<I>(rhs))...
            );
        }

        template <typename MatrixLHS, typename MatrixRHS, std::size_t ...I>
        auto element_diff_impl (MatrixLHS lhs,
                                MatrixRHS rhs,
                                std::index_sequence<I...>)
        {
            return hana::make<hana::Tuple>(
                (tuple_access::get<I>(lhs) - tuple_access::get<I>(rhs))...
            );
        }

        template <typename MatrixLHS, typename MatrixRHS, std::size_t ...I>
        auto element_prod_impl (MatrixLHS lhs,
                                MatrixRHS rhs,
                                std::index_sequence<I...>)
        {
            return hana::make<hana::Tuple>(
                (tuple_access::get<I>(lhs) * tuple_access::get<I>(rhs))...
            );
        }

        template <typename MatrixLHS, typename MatrixRHS, std::size_t ...I>
        auto element_div_impl (MatrixLHS lhs,
                               MatrixRHS rhs,
                               std::index_sequence<I...>)
        {
            return hana::make<hana::Tuple>(
                (tuple_access::get<I>(lhs) / tuple_access::get<I>(rhs))...
            );
        }


#if 0
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
#endif

        template <typename Matrix, std::size_t ...I>
        auto neg_impl (Matrix m, std::index_sequence<I...>)
        { return hana::make<hana::Tuple>(-tuple_access::get<I>(m)...); }

        template <typename Tuple, bool Abs>
        struct sum_impl
        {
            template <std::size_t I, typename T, typename State>
            auto call (State prev)
            {
                auto value = hana::at_c<I>(t_);
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
                auto value = hana::at_c<I>(t_);
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
                using common_type = decltype(prev.second + hana::at_c<I>(t_));
                using retval_type = std::pair<std::size_t, common_type>;
                using std::abs;
                auto value = abs(hana::at_c<I>(t_));
                return
                    prev.second < value ?
                    retval_type{I, value} :
                    retval_type{prev.first, prev.second};
            }

            Tuple t_;
        };


        // matrix determinant type
#if BOOST_UNITS_BLAS_USE_INEXACT_DETERMINANT_TYPE
        template <typename Matrix, std::size_t I, std::size_t N>
        struct simplified_determinant_type_impl
        {
            template <typename Prev>
            static constexpr auto call (Prev prev)
            {
                using type = tuple_element_t<
                    I * Matrix::num_columns + I,
                    Matrix
                >;
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
            using first = tuple_element_t<0, Matrix>;
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
        struct determinant_type;

        template <typename Matrix>
        using determinant_type_t = typename determinant_type<Matrix>::type;

        template <std::size_t N, typename InSeq, typename OutSeq>
        struct filter_out_n;

        template <std::size_t N,
                  typename OutSeq,
                  std::size_t Head,
                  std::size_t ...Tail>
        struct filter_out_n<N, std::index_sequence<Head, Tail...>, OutSeq>
        {
            using type = typename filter_out_n<
                N,
                std::index_sequence<Tail...>,
                std::conditional_t<
                    Head == N,
                    OutSeq,
                    decltype(push_back<Head>(OutSeq{}))
                >
            >::type;
        };

        template <std::size_t N, typename OutSeq>
        struct filter_out_n<N, std::index_sequence<>, OutSeq>
        {
            using type = OutSeq;
        };

        template <typename Rows, std::size_t Column, typename Matrix>
        constexpr auto subdeterminant (Matrix m)
        {
            using columns = typename filter_out_n<
                Column,
                decltype(std::make_index_sequence<Matrix::num_columns>()),
                std::index_sequence<>
            >::type;
            constexpr auto seqs = slice_indices_and_types<Matrix, columns>(
                std::pair<std::index_sequence<>, type_sequence<>>{},
                Rows{}
            );
            using minor_matrix = matrix_t<
                decltype(tuple_from_types(seqs.second)),
                Matrix::num_rows - 1,
                Matrix::num_columns - 1
            >;
            return determinant_type_t<minor_matrix>{};
        }

        struct sum_type
        {
            template <typename T1, typename T2>
            struct apply
            {
                using type = decltype(std::declval<T1>() + std::declval<T2>());
            };
        };

        template <typename Matrix, std::size_t Head, std::size_t ...Tail>
        constexpr auto
        determinant_type_impl (Matrix m, std::index_sequence<Head, Tail...>)
        {
            using rows = std::index_sequence<Tail...>;
            using tuple = hana::_tuple<
                decltype(tuple_access::get<Head>(m) *
                         subdeterminant<rows, Head>(m)),
                decltype(tuple_access::get<Tail>(m) *
                         subdeterminant<rows, Tail>(m))...
            >;

            using sum = static_foldl_t<
                sum_type,
                typename std::tuple_element<0, tuple>::type,
                decltype(types_from_tuple(tuple{}))
            >;

            return sum{};
        }

        template <typename Matrix>
        struct determinant_type
        {
            using type = decltype(
                determinant_type_impl(
                    Matrix{},
                    std::make_index_sequence<Matrix::num_columns>()
                )
            );
        };

        template <typename Tuple>
        struct determinant_type<matrix_t<Tuple, 1, 1>>
        {
            using type = decltype(determinant(matrix_t<Tuple, 1, 1>{}));
        };

        template <typename Tuple>
        struct determinant_type<matrix_t<Tuple, 2, 2>>
        {
            using type = decltype(determinant(matrix_t<Tuple, 2, 2>{}));
        };

        template <typename Tuple>
        struct determinant_type<matrix_t<Tuple, 3, 3>>
        {
            using type = decltype(determinant(matrix_t<Tuple, 3, 3>{}));
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
                    tuple_element_t<I, Matrix>,
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
                    tuple_element_t<I, Matrix>,
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
                    tuple_element_t<I, Matrix>,
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
                    tuple_element_t<I, Matrix>,
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
        matrix<>.  Both of @c RowIndices and @c ColumnIndices must be a @c
        std::index_sequence<> ; all values in @c RowIndices and @c ColumnIndices
        must be less than the number of rows and columns in @c m, respectively.
        Note that duplication and order preservation are not enforced for the
        values in @c RowIndices and @c ColumnIndices.  It is therefore possible
        to use @c slice<>() to rearrange and/or duplicate rows and/or
        columns. */
    template <typename RowIndices,
              typename ColumnIndices,
              typename Tuple,
              std::size_t Rows,
              std::size_t Columns>
    auto slice (matrix_t<Tuple, Rows, Columns> m)
    {
        using matrix_type = matrix_t<Tuple, Rows, Columns>;

        static_assert(
            0 < RowIndices::size(),
            "slice() requires at least one row index to compute its result"
        );

        static_assert(
            0 < ColumnIndices::size(),
            "slice() requires at least one column index to compute its result"
        );

        auto seqs = detail::slice_indices_and_types<matrix_type, ColumnIndices>(
            std::pair<std::index_sequence<>, detail::type_sequence<>>{},
            RowIndices{}
        );

        auto tuple = detail::tuple_from_types(seqs.second);
        iterate_indexed(
            detail::tuple_assign<decltype(tuple), matrix_type>{tuple, m},
            seqs.first
        );

        return detail::make_matrix<RowIndices::size(), ColumnIndices::size()>(
            tuple
        );
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
        auto seqs = detail::transpose_indices_and_types<matrix_type>(
            std::make_index_sequence<matrix_type::num_elements>()
        );
        auto tuple = tuple_from_types(seqs.second);
        iterate_indexed(
            detail::tuple_assign<decltype(tuple), matrix_type>{tuple, m},
            seqs.first
        );
        return detail::make_matrix<Columns, Rows>(tuple);
    }

    /** Returns the negation of @c m.  @c m must be a @c matrix<>. */
    template <std::size_t Rows, std::size_t Columns, typename ...T>
    auto neg (matrix_t<hana::_tuple<T...>, Rows, Columns> m)
    {
        return detail::make_matrix<Rows, Columns>(
            detail::neg_impl(
                m,
                std::make_index_sequence<Rows * Columns>()
            )
        );
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
    {
        return detail::make_matrix<Rows, Columns>(
            detail::element_sum_impl(
                lhs,
                rhs,
                std::make_index_sequence<Rows * Columns>()
            )
        );
    }

    /** Returns the elementwise difference of @c lhs and @c rhs. @c lhs and @c
        rhs must be <c>matrix<></c>s with the same dimensions.  Also, every
        difference <c>lhs(i, j) - rhs(i, j)</c> must be a valid operation.  */
    template <typename Tuple1,
              typename Tuple2,
              std::size_t Rows,
              std::size_t Columns>
    auto diff (matrix_t<Tuple1, Rows, Columns> lhs,
               matrix_t<Tuple2, Rows, Columns> rhs)
    {
        return detail::make_matrix<Rows, Columns>(
            detail::element_diff_impl(
                lhs,
                rhs,
                std::make_index_sequence<Rows * Columns>()
            )
        );
    }

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
        return detail::make_matrix<Rows1, Columns2>(
            detail::matrix_prod_impl(
                lhs,
                rhs,
                std::make_index_sequence<Rows1 * Columns2>()
            )
        );
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
        return detail::make_matrix<Rows, Columns>(
            detail::scalar_prod_impl(
                m,
                t,
                std::make_index_sequence<Rows * Columns>()
            )
        );
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
        return detail::make_matrix<Rows, Columns>(
            detail::scalar_div_impl(
                m,
                t,
                std::make_index_sequence<Rows * Columns>()
            )
        );
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
        return detail::make_matrix<Rows, Columns>(
            detail::element_prod_impl(
                lhs,
                rhs,
                std::make_index_sequence<Rows * Columns>()
            )
        );
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
        return detail::make_matrix<Rows, Columns>(
            detail::element_div_impl(
                lhs,
                rhs,
                std::make_index_sequence<Rows * Columns>()
            )
        );
    }

    /** Swaps the values in @c lhs and @c rhs.  Note that this is an
        O(@c size<matrix<T> >::value) operation. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    void swap (matrix_t<Tuple, Rows, Columns> & lhs,
               matrix_t<Tuple, Rows, Columns> & rhs)
    {
        using std::swap;
        hana::fold.left(
            detail::tuple_access::all(lhs),
            hana::size_t<0>,
            [&](auto i, auto & x) {
                swap(hana::at(detail::tuple_access::get<i.value>(rhs), i), x);
                return hana::succ(i);
            }
        );
#if 0
        using matrix_type = matrix_t<Tuple, Rows, Columns>;
        detail::iterate_simple<matrix_type::num_elements>(
            detail::swap<matrix_type>{lhs, rhs}
        );
#endif
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
            hana::at_c<1>(l) * hana::at_c<2>(r) -
            hana::at_c<2>(l) * hana::at_c<1>(r);
        auto _1 =
            hana::at_c<2>(l) * hana::at_c<0>(r) -
            hana::at_c<0>(l) * hana::at_c<2>(r);
        auto _2 =
            hana::at_c<0>(l) * hana::at_c<1>(r) -
            hana::at_c<1>(l) * hana::at_c<0>(r);
        return detail::make_matrix<Rows, Columns>(hana::make<hana::Tuple>(_0, _1, _2));
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
    auto sum (matrix_t<hana::_tuple<Head, Tail...>, Rows, Columns> v,
              typename std::enable_if<
                  Rows == 1 || Columns == 1
              >::type* = 0)
    {
        auto f = detail::sum_impl<hana::_tuple<Head, Tail...>, false>{
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
    auto norm_1 (matrix_t<hana::_tuple<Head, Tail...>, Rows, Columns> v,
                 typename std::enable_if<
                     Rows == 1 || Columns == 1
                 >::type* = 0)
    {
        auto f = detail::sum_impl<hana::_tuple<Head, Tail...>, true>{
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
    auto norm_2 (matrix_t<hana::_tuple<Head, Tail...>, Rows, Columns> v,
                 typename std::enable_if<
                     Rows == 1 || Columns == 1
                 >::type* = 0)
    {
        auto f = detail::norm_2_impl<hana::_tuple<Head, Tail...>>{
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
    auto norm_inf (matrix_t<hana::_tuple<T...>, Rows, Columns> v,
                   typename std::enable_if<
                       Rows == 1 || Columns == 1
                   >::type* = 0)
    {
        auto f = detail::norm_inf_index_impl<hana::_tuple<T...>>{
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
    auto norm_inf_index (matrix_t<hana::_tuple<T...>, Rows, Columns> v,
                         typename std::enable_if<
                             Rows == 1 || Columns == 1
                         >::type* = 0)
    {
        auto f = detail::norm_inf_index_impl<hana::_tuple<T...>>{
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

        using result_type = detail::determinant_type_t<matrix_type>;
        using raw_value_type = detail::value_type_t<result_type>;
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

        using result_type = detail::inverse_type_t<matrix_type>;
        using temp_value_type = detail::determinant_type_t<matrix_type>;
        using raw_value_type = detail::value_type_t<temp_value_type>;
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

        using temp_value_type = detail::determinant_type_t<a_type>;
        using raw_value_type = detail::value_type_t<temp_value_type>;

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
