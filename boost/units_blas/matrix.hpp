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

        struct tuple_access
        {
            template <std::size_t I, typename Matrix>
            static decltype(auto) get (Matrix && m)
            { return std::get<I>(m.data_); }

            template <typename Matrix>
            static decltype(auto) all (Matrix && m)
            { return (m.data_); }
        };

        template <std::size_t I, std::size_t N, typename F>
        struct iterate_simple_impl;

        template <std::size_t I, std::size_t N, typename F>
        struct iterate_simple_impl
        {
            static void call (F f)
            {
                f.template call<I>();
                iterate_simple_impl<I + 1, N, F>::call(f);
            }
        };

        template <std::size_t N, typename F>
        struct iterate_simple_impl<N, N, F>
        {
            static void call (F f) {}
        };

        template <std::size_t N, typename F>
        void iterate_simple (F f)
        { iterate_simple_impl<0, N, F>::call(f); }

        template <typename MatrixL, typename MatrixR>
        struct assign
        {
            template <std::size_t I>
            void call ()
            { tuple_access::get<I>(lhs) = tuple_access::get<I>(rhs); }

            MatrixL & lhs;
            MatrixR const & rhs;
        };

        template <typename MatrixL, typename MatrixR>
        struct plus_assign
        {
            template <std::size_t I>
            void call ()
            { tuple_access::get<I>(lhs) += tuple_access::get<I>(rhs); }

            MatrixL & lhs;
            MatrixR const & rhs;
        };

        template <typename MatrixL, typename MatrixR>
        struct minus_assign
        {
            template <std::size_t I>
            void call ()
            { tuple_access::get<I>(lhs) -= tuple_access::get<I>(rhs); }

            MatrixL & lhs;
            MatrixR const & rhs;
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
        matrix_t () :
            data_ ()
        {}

        /** Copy ctor. */
        matrix_t (matrix_t const & m) :
            data_ (m.data_)
        {}

        /** Copy ctor.  @c *this and @c m must have the same dimensions, and
            every @c m(i, j) must be convertible to the type of <c>(*this)(i,
            j)</c>. */
        template <typename MatrixR>
        matrix_t (MatrixR const & m,
                typename std::enable_if<
                    MatrixR::num_rows == num_rows &&
                    MatrixR::num_columns == num_columns
                >::type * = 0) :
            data_ ()
        { assign_data(m); }

        /** Assignment operator. */
        matrix_t & operator= (matrix_t const & rhs)
        {
            data_ = rhs.data_;
            return *this;
        }

        /** Assignment operator.  @c *this and @c m must have the same
            dimensions, and every @c m(i, j) must be convertible to the type
            of <c>(*this)(i, j)</c>. */
        template <typename MatrixR>
        typename std::enable_if<
            MatrixR::num_rows == num_rows &&
            MatrixR::num_columns == num_columns,
            matrix_t &
        >::type
        operator= (MatrixR const & rhs)
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

#if 0
        /** Prints @c *this to stream @c os.  Note that this function remains
            unimplemented unless @c boost/units_blas/io.hpp is also
            included. */
        std::ostream & print (std::ostream & os) const;
#endif


        /** Returns a non-const reference to the element at row @c I, column @c
            J. */
        template <size_type I, size_type J>
        decltype(auto) at ()
        { return std::get<I * num_columns + J>(data_); }

        /** Adds the values in @c rhs to @c *this, element-by-element.  @c
            *this and @c m must have the same dimensions, and every sum
            <c>(*this)(i, j) + rhs(i, j)</c> must be convertible to the type
            of <c>(*this)(i, j)</c>. */
        template <typename MatrixR>
        typename std::enable_if<
            MatrixR::num_rows == num_rows &&
            MatrixR::num_columns == num_columns,
            matrix_t &
        >::type
        operator+= (MatrixR const & rhs)
        {
            detail::iterate_simple<num_elements>(
                detail::plus_assign<self_type, MatrixR>{*this, rhs}
            );
            return *this;
        }

        /** Subtracts the values in @c rhs from @c *this,
            element-by-element. @c *this and @c m must have the same
            dimensions, and every difference <c>(*this)(i, j) - rhs(i, j)</c>
            must be convertible to the type of <c>(*this)(i, j)</c>. */
        template <typename MatrixR>
        typename std::enable_if<
            MatrixR::num_rows == num_rows &&
            MatrixR::num_columns == num_columns,
            matrix_t &
        >::type
        operator-= (MatrixR const & rhs)
        {
            detail::iterate_simple<num_elements>(
                detail::minus_assign<self_type, MatrixR>{*this, rhs}
            );
            return *this;
        }

        /** Multiplies each element of @c *this by @c val.  Every product
            <c>(*this)(i, j) * val</c> must be convertible to the type of
            <c>(*this)(i, j)</c>. */
        template <typename T>
        matrix_t & operator*= (T const & val)
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
        matrix_t & operator/= (T const & val)
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

#ifndef BOOST_UNITS_BLAS_DOXYGEN
#if 0 // TODO
        struct column_lengths_equal :
            fusion::result_of::equal_to<
                typename fusion::result_of::find_if<
                    value_types,
                    mpl::not_equal_to<
                        fusion::result_of::size<mpl::_1>,
                        mpl::int_<num_columns_t::value>
                    >
                >::type,
                typename fusion::result_of::end<value_types>::type
            >
        {};

        // If you're seeing an error here, it's because you defined two or
        // more rows that have different numbers of columns.
        BOOST_MPL_ASSERT((column_lengths_equal));
#endif

        static_assert(
            std::tuple_size<Tuple>::value == num_elements,
            "matrix_t<> must be specified with a Tuple template parameter "
            "consisting of Rows * Columns types"
        );
        static_assert(
            0 < num_rows,
            "matrix_t<> must be specified with a positive Rows template parameter"
        );
        static_assert(
            0 < num_columns,
            "matrix_t<> must be specified with a positive Columns template parameter"
        );

        template <typename MatrixR>
        void assign_data (MatrixR const & rhs)
        {
            detail::iterate_simple<num_elements>(
                detail::assign<self_type, MatrixR>{*this, rhs}
            );
        }
#endif

        value_types data_;
    };

    namespace detail {

        template <typename HeadRow, typename ...TailRows>
        struct matrix_type
        {
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
