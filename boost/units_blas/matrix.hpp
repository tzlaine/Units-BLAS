// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_MATRIX_HPP
#define BOOST_UNITS_BLAS_MATRIX_HPP

#include <boost/units_blas/config.hpp>
#include <boost/units_blas/matrix_fwd.hpp>

#include <tuple>

#include <iosfwd>


namespace boost { namespace units_blas {

    namespace detail {

        // tuple access
        struct tuple_access
        {
            template <std::size_t I, typename Matrix>
            static decltype(auto) get (Matrix && m)
            { return std::get<I>(m.data_); }

            template <typename Matrix>
            static decltype(auto) all (Matrix && m)
            { return (m.data_); }
        };


        // sequences
        template <typename ...T>
        struct type_sequence
        {
            static const std::size_t size = sizeof...(T);
        };


        // push_back
        template <std::size_t Tail, std::size_t ...Head>
        constexpr auto push_back (std::index_sequence<Head...>)
        { return std::index_sequence<Head..., Tail>{}; }

        template <typename Tail, typename ...Head>
        constexpr auto push_back (type_sequence<Head...>)
        { return type_sequence<Head..., Tail>{}; }


        // tuple <--> typelist

        // This is not declared as return-auto as a workaround for a bug in
        // Clang 3.4.
        template <typename ...T>
        std::tuple<T...> tuple_from_types (type_sequence<T...>)
        { return std::tuple<T...>{}; }

        template <typename ...T>
        constexpr auto types_from_tuple (std::tuple<T...>)
        { return type_sequence<T...>{}; }


        // matrix from tuple
        template <std::size_t Rows, std::size_t Columns, typename Tuple>
        auto make_matrix (Tuple t)
        {
            matrix_t<Tuple, Rows, Columns> retval;
            tuple_access::all(retval) = t;
            return retval;
        }


        template <typename Matrix>
        constexpr std::size_t transpose_index (std::size_t i)
        {
            return
                i % Matrix::num_columns * Matrix::num_rows +
                i / Matrix::num_columns;
        }


        // tuple element I of a matrix_t
        template <std::size_t I, typename Matrix>
        struct tuple_element :
            std::tuple_element<I, typename Matrix::value_types>
        {};

        template <std::size_t I, typename Matrix>
        using tuple_element_t = typename tuple_element<I, Matrix>::type;


        // iterate simple
        void iterate_landing_pad (...) {}

        template <std::size_t I, typename F>
        int iterate_fn_wrapper (F f)
        {
            f.template call<I>();
            return 0;
        }

        template <typename F, std::size_t ...I>
        void iterate_simple_impl (F f, std::index_sequence<I...>)
        { iterate_landing_pad(iterate_fn_wrapper<I>(f)...); }

        template <std::size_t N, typename F>
        void iterate_simple (F f)
        { iterate_simple_impl(f, std::make_index_sequence<N>()); }


        // iterate_simple function objects

        template <typename MatrixL, typename MatrixR>
        struct assign
        {
            template <std::size_t I>
            void call ()
            { tuple_access::get<I>(lhs) = tuple_access::get<I>(rhs); }

            MatrixL & lhs;
            MatrixR rhs;
        };

        template <typename MatrixL, typename MatrixR>
        struct plus_assign
        {
            template <std::size_t I>
            void call ()
            { tuple_access::get<I>(lhs) += tuple_access::get<I>(rhs); }

            MatrixL & lhs;
            MatrixR rhs;
        };

        template <typename MatrixL, typename MatrixR>
        struct minus_assign
        {
            template <std::size_t I>
            void call ()
            { tuple_access::get<I>(lhs) -= tuple_access::get<I>(rhs); }

            MatrixL & lhs;
            MatrixR rhs;
        };

        template <typename Matrix, typename T>
        struct scalar_multiply_assign
        {
            template <std::size_t I>
            void call ()
            { tuple_access::get<I>(lhs) *= val; }

            Matrix & lhs;
            T val;
        };

        template <typename Matrix, typename T>
        struct scalar_divide_assign
        {
            template <std::size_t I>
            void call ()
            { tuple_access::get<I>(lhs) /= val; }

            Matrix & lhs;
            T val;
        };

    }

    /** Note that in the member documentation that follows, the notation @c
        m(i, j) is used to refer to the element of matrix @c m at row @c i,
        column @c j, even though there is no @c operator() defined for @c
        matrix_t<>. */
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    class matrix_t
    {
    public:
        using self_type = matrix_t<Tuple, Rows, Columns>;
        using size_type = std::size_t;
        using value_types = Tuple;

        static const std::size_t num_rows = Rows;
        static const std::size_t num_columns = Columns;
        static const std::size_t num_elements = num_rows * num_columns;

        /** Default ctor.  Each element of @c *this is default constructed. */
        matrix_t () = default;

        /** Copy ctor. */
        matrix_t (matrix_t const & m) = default;

        /** Copy ctor.  @c *this and @c m must have the same dimensions, and
            every @c m(i, j) must be convertible to the type of <c>(*this)(i,
            j)</c>. */
        template <typename Tuple2>
        matrix_t (matrix_t<Tuple2, Rows, Columns> m) :
            data_ ()
        { assign_data(m); }

        /** Assignment operator. */
        matrix_t & operator= (matrix_t const & rhs) = default;

        /** Assignment operator.  @c *this and @c m must have the same
            dimensions, and every @c m(i, j) must be convertible to the type
            of <c>(*this)(i, j)</c>. */
        template <typename TupleRHS>
        matrix_t & operator= (matrix_t<TupleRHS, Rows, Columns> rhs)
        {
            assign_data(rhs);
            return *this;
        }

        /** Returns the number of elements in @c *this. */
        size_type size () const
        { return num_elements; }

        /** Returns the number of rows in @c *this. */
        size_type rows () const
        { return num_rows; }

        /** Returns the number of columns in @c *this. */
        size_type columns () const
        { return num_columns; }

        /** Returns a const reference to the element at row @c I, column @c
            J. */
        template <size_type I, size_type J>
        auto at () const
        { return std::get<I * num_columns + J>(data_); }

        /** Prints @c *this to stream @c os.  Note that this function remains
            unimplemented unless @c boost/units_blas/io.hpp is also
            included. */
        std::ostream & print (std::ostream & os) const;

        /** Returns a non-const reference to the element at row @c I, column @c
            J. */
        template <size_type I, size_type J>
        decltype(auto) at ()
        { return std::get<I * num_columns + J>(data_); }

        /** Adds the values in @c rhs to @c *this, element-by-element.  @c
            *this and @c m must have the same dimensions, and every sum
            <c>(*this)(i, j) + rhs(i, j)</c> must be convertible to the type
            of <c>(*this)(i, j)</c>. */
        template <typename TupleRHS>
        matrix_t & operator+= (matrix_t<TupleRHS, Rows, Columns> rhs)
        {
            using matrix_rhs = matrix_t<TupleRHS, Rows, Columns>;
            detail::iterate_simple<num_elements>(
                detail::plus_assign<self_type, matrix_rhs>{*this, rhs}
            );
            return *this;
        }

        /** Subtracts the values in @c rhs from @c *this,
            element-by-element. @c *this and @c m must have the same
            dimensions, and every difference <c>(*this)(i, j) - rhs(i, j)</c>
            must be convertible to the type of <c>(*this)(i, j)</c>. */
        template <typename TupleRHS>
        matrix_t & operator-= (matrix_t<TupleRHS, Rows, Columns> rhs)
        {
            using matrix_rhs = matrix_t<TupleRHS, Rows, Columns>;
            detail::iterate_simple<num_elements>(
                detail::minus_assign<self_type, matrix_rhs>{*this, rhs}
            );
            return *this;
        }

        /** Multiplies each element of @c *this by @c val.  Every product
            <c>(*this)(i, j) * val</c> must be convertible to the type of
            <c>(*this)(i, j)</c>. */
        template <typename T>
        matrix_t & operator*= (T val)
        {
            detail::iterate_simple<num_elements>(
                detail::scalar_multiply_assign<self_type, T>{*this, val}
            );
            return *this;
        }

        /** Divides each element of @c *this by @c val.  Every quotient
            <c>(*this)(i, j) / val</c> must be convertible to the type of
            <c>(*this)(i, j)</c>. */
        template <typename T>
        matrix_t & operator/= (T val)
        {
            detail::iterate_simple<num_elements>(
                detail::scalar_divide_assign<self_type, T>{*this, val}
            );
            return *this;
        }

    private:
        template <typename Tuple_, std::size_t Rows_, std::size_t Columns_>
        friend class matrix_t;

        friend struct detail::tuple_access;

        static_assert(
            std::tuple_size<Tuple>::value == num_elements,
            "matrix_t<> must be specified with a Tuple template parameter "
            "consisting of Rows * Columns types"
        );

        static_assert(
            0 < num_rows,
            "matrix_t<> must be specified with a positive Rows "
            "template parameter"
        );

        static_assert(
            0 < num_columns,
            "matrix_t<> must be specified with a positive Columns "
            "template parameter"
        );

        template <typename MatrixR>
        void assign_data (MatrixR rhs)
        {
            detail::iterate_simple<num_elements>(
                detail::assign<self_type, MatrixR>{*this, rhs}
            );
        }

        value_types data_;
    };

    namespace detail {

        template <typename HeadRow, typename ...TailRows>
        constexpr auto tuples_same_size (bool same,
                                         std::size_t prev,
                                         type_sequence<HeadRow, TailRows...>)
        {
            constexpr std::size_t size = std::tuple_size<HeadRow>::value;
            return tuples_same_size(
                same && prev == size,
                size,
                type_sequence<TailRows...>{}
            );
        }

        constexpr auto tuples_same_size (bool s, std::size_t, type_sequence<>)
        { return s; }

        template <typename HeadRow, typename ...TailRows>
        struct matrix_type
        {
            static_assert(
                tuples_same_size(
                    true,
                    std::tuple_size<HeadRow>::value,
                    type_sequence<HeadRow, TailRows...>{}
                ),
                "matrix_type<> requires tuples of uniform length"
            );

            using tuple = decltype(std::tuple_cat(
                std::declval<HeadRow>(),
                std::declval<TailRows>()...
            ));
            static const std::size_t rows = sizeof...(TailRows) + 1;
            static const std::size_t columns = std::tuple_size<HeadRow>::value;
            using type = matrix_t<tuple, rows, columns>;
        };

    }

    template <typename ...Rows>
    using matrix = typename detail::matrix_type<Rows...>::type;

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MATRIX_HPP
