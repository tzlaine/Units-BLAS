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

#include <boost/hana.hpp>

#include <tuple>

#include <iosfwd>


namespace boost { namespace units_blas {

    namespace detail {

        // tuple access
        struct tuple_access
        {
            template <std::size_t I, typename Matrix>
            static decltype(auto) get (Matrix && m)
            { return hana::at_c<I>(m.data_); }

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
        hana::_tuple<T...> tuple_from_types (type_sequence<T...>)
        { return hana::_tuple<T...>{}; }

        template <typename ...T>
        constexpr auto types_from_tuple (hana::_tuple<T...>)
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
                i % Matrix::num_rows * Matrix::num_columns +
                i / Matrix::num_rows;
        }


        // tuple element I of a matrix_t
        template <std::size_t I, typename Matrix>
        using tuple_element_t =
            std::remove_reference_t<
                decltype(
                    hana::at_c<I>(std::declval<typename Matrix::value_types>())
                )
            >;


#if 1
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
#endif

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
        matrix_t () {}

        /** Copy ctor. */
        matrix_t (matrix_t const & m) :
            data_ (m.data_)
        {}

        /** Copy ctor.  @c *this and @c m must have the same dimensions, and
            every @c m(i, j) must be convertible to the type of <c>(*this)(i,
            j)</c>. */
        template <typename Tuple2>
        matrix_t (matrix_t<Tuple2, Rows, Columns> m) :
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
        { return hana::at_c<I * num_columns + J>(data_); }

        /** Prints @c *this to stream @c os.  Note that this function remains
            unimplemented unless @c boost/units_blas/io.hpp is also
            included. */
        std::ostream & print (std::ostream & os) const;

        /** Returns a non-const reference to the element at row @c I, column @c
            J. */
        template <size_type I, size_type J>
        decltype(auto) at ()
        { return hana::at_c<I * num_columns + J>(data_); }

        /** Adds the values in @c rhs to @c *this, element-by-element.  @c
            *this and @c m must have the same dimensions, and every sum
            <c>(*this)(i, j) + rhs(i, j)</c> must be convertible to the type
            of <c>(*this)(i, j)</c>. */
        template <typename TupleRHS>
        matrix_t & operator+= (matrix_t<TupleRHS, Rows, Columns> rhs)
        {
            hana::fold.left(
                rhs.data_,
                hana::size_t<0>,
                [&](auto i, auto x) {
                    hana::at(data_, i) += x;
                    return hana::succ(i);
                }
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
            hana::fold.left(
                rhs.data_,
                hana::size_t<0>,
                [&](auto i, auto x) {
                    hana::at(data_, i) -= x;
                    return hana::succ(i);
                }
            );
            return *this;
        }

        /** Multiplies each element of @c *this by @c val.  Every product
            <c>(*this)(i, j) * val</c> must be convertible to the type of
            <c>(*this)(i, j)</c>. */
        template <typename T>
        matrix_t & operator*= (T val)
        {
            hana::fold.left(
                data_,
                hana::size_t<0>,
                [&](auto i, auto) {
                    hana::at(data_, i) *= val;
                    return hana::succ(i);
                }
            );
            return *this;
        }

        /** Divides each element of @c *this by @c val.  Every quotient
            <c>(*this)(i, j) / val</c> must be convertible to the type of
            <c>(*this)(i, j)</c>. */
        template <typename T>
        matrix_t & operator/= (T val)
        {
            hana::fold.left(
                data_,
                hana::size_t<0>,
                [&](auto i, auto) {
                    hana::at(data_, i) /= val;
                    return hana::succ(i);
                }
            );
            return *this;
        }

    private:
        template <typename Tuple_, std::size_t Rows_, std::size_t Columns_>
        friend class matrix_t;

        friend struct detail::tuple_access;

        static_assert(
            Tuple::size == num_elements,
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
#if 0

            // Could not get this to work due to 'static_assert failed
            // "hana::zip.shortest.with(f, xs, ys...) requires ys to be a
            // Sequence'

            // Awkward (and maybe inefficient?) to use zip.with for this;
            // would prefer to use some kind of mutating copy.
            data_ = hana::zip.with(rhs.data_, data_, [](auto x, auto y) {
                return static_cast<decltype(y)>(x);
            });
#elif 0
            // Close, but also does not work, as the sequence is not Foldable.
            hana::for_each(
                hana::detail::generate_index_sequence<Rows * Columns>{},
                [&](auto i) {
                    hana::at(data_, i) =
                        static_cast<std::remove_reference_t<decltype(hana::at(data_, i))>>(
                            hana::at(rhs.data_, i)
                        );
                }
            );
#elif 1
            // Feels dirty, but works
            hana::fold.left(
                rhs.data_,
                hana::size_t<0>,
                [&](auto i, auto x) {
                    hana::at(data_, i) =
                        static_cast<std::remove_reference_t<decltype(hana::at(data_, i))>>(x);
                    return hana::succ(i);
                }
            );
#endif
        }

        value_types data_;
    };

    namespace detail {

        template <typename HeadRow, typename ...TailRows>
        struct matrix_type
        {
            static const std::size_t rows = sizeof...(TailRows) + 1;
#if 0 // attempt 1
            static const std::size_t columns = hana::size(HeadRow{});
#elif 0 // attempt 2
            static const std::size_t columns = hana::size(std::declval<HeadRow>());
#elif 0 // attempt 3
            static const std::size_t columns = hana::size(hana::type<HeadRow>);
#else // This compiles.  Sigh.
            static const std::size_t columns = HeadRow::size;
#endif

            struct same_size_as_head
            {
                template <typename T>
                constexpr bool operator() (hana::_type<T> x)
                { return hana::_type<T>::_::type::size == columns;}
            };

            static_assert(
                hana::all_of(
                    hana::make<hana::Tuple>(hana::_type<TailRows>{}...),
                    same_size_as_head{}
                ),
                "matrix_type<> requires tuples of uniform length"
            );

            using tuple = decltype(hana::flatten(
                hana::make<hana::Tuple>(
                    HeadRow{},
                    TailRows{}...
                )
            ));

            using type = matrix_t<tuple, rows, columns>;
        };

    }

    template <typename ...Rows>
    using matrix = typename detail::matrix_type<Rows...>::type;

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MATRIX_HPP
