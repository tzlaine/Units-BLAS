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
#include <boost/units_blas/iterate.hpp>
#include <boost/units_blas/result_of/at.hpp>
#include <boost/units_blas/detail/simple_iteration.hpp>

#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/not_equal_to.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/mpl/times.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/fusion/algorithm/query/find_if.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/mpl.hpp>
#include <boost/fusion/sequence/comparison/equal_to.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/fusion/sequence/intrinsic/value_at.hpp>

#include <iosfwd>


namespace boost { namespace units_blas {

    namespace detail {
        template <typename Rows>
        struct deep_as_vector
        {
            typedef typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    Rows,
                    fusion::result_of::as_vector<mpl::_1>
                >::type
            >::type type;
        };

        template <typename ValueTypes>
        struct num_rows
        {
            typedef mpl::size_t<
                fusion::result_of::size<ValueTypes>::type::value
            > type;
        };

        template <typename ValueTypes>
        struct num_columns
        {
            typedef mpl::size_t<
                fusion::result_of::size<
                    typename fusion::result_of::value_at_c<
                        ValueTypes,
                        0
                    >::type
                >::type::value
            > type;
        };
    }

    /** Note that in the member documentation that follows, the notation @c
        m(i, j) is used to refer to the element of matrix @c m at row @c i,
        column @c j, even though there is no @c operator() defined for @c
        matrix<>. */
    template <typename Rows>
    class matrix
    {
    public:
        typedef std::size_t size_type;
        typedef typename detail::deep_as_vector<Rows>::type value_types;
        typedef typename detail::num_rows<value_types>::type num_rows_t;
        typedef typename detail::num_columns<value_types>::type num_columns_t;

        /** Default ctor.  Each element of @c *this is default constructed. */
        matrix () :
            data_ ()
            {}

        /** Copy ctor. */
        matrix (matrix const & m) :
            data_ (m.data_)
            {}

        /** Copy ctor.  @c *this and @c m must have the same dimensions, and
            every @c m(i, j) must be convertible to the type of <c>(*this)(i,
            j)</c>. */
        template <typename T>
        matrix (matrix<T> const & m,
                typename enable_if<
                    mpl::and_<
                        mpl::equal_to<typename matrix<T>::num_rows_t, num_rows_t>,
                        mpl::equal_to<typename matrix<T>::num_columns_t, num_columns_t>
                    >
                >::type * = 0) :
            data_ ()
            { assign_data(m); }

        /** Assignment operator. */
        matrix & operator= (matrix const & rhs)
            {
                data_ = rhs.data_;
                return *this;
            }

        /** Assignment operator.  @c *this and @c m must have the same
            dimensions, and every @c m(i, j) must be convertible to the type
            of <c>(*this)(i, j)</c>. */
        template <typename T>
        typename enable_if<
            mpl::and_<
                mpl::equal_to<typename matrix<T>::num_rows_t, num_rows_t>,
                mpl::equal_to<typename matrix<T>::num_columns_t, num_columns_t>
            >,
            matrix &
        >::type
        operator= (matrix<T> const & rhs)
            {
                assign_data(rhs);
                return *this;
            }

        /** Returns the number of elements in @c *this. */
        size_type size () const
            { return mpl::times<num_rows_t, num_columns_t>::value; }

        /** Returns the number of rows in @c *this. */
        size_type rows () const
            { return num_rows_t::value; }

        /** Returns the number of columns in @c *this. */
        size_type columns () const
            { return num_columns_t::value; }

        /** Returns a const reference to the element at row @c I, column @c
            J. */
        template <size_type I, size_type J>
        typename result_of::at_c<matrix const, I, J>::type
        at () const
            { return fusion::at_c<J>(fusion::at_c<I>(data_)); }

        /** Prints @c *this to stream @c os.  Note that this function remains
            unimplemented unless @c boost/units_blas/io.hpp is also
            included. */
        std::ostream & print (std::ostream & os) const;


        /** Returns a non-const reference to the element at row @c I, column @c
            J. */
        template <size_type I, size_type J>
        typename result_of::at_c<matrix, I, J>::type
        at ()
            { return fusion::at_c<J>(fusion::at_c<I>(data_)); }

        /** Adds the values in @c rhs to @c *this, element-by-element.  @c
            *this and @c m must have the same dimensions, and every sum
            <c>(*this)(i, j) + rhs(i, j)</c> must be convertible to the type
            of <c>(*this)(i, j)</c>. */
        template <typename T>
        typename enable_if<
            mpl::and_<
                mpl::equal_to<typename matrix<T>::num_rows_t, num_rows_t>,
                mpl::equal_to<typename matrix<T>::num_columns_t, num_columns_t>
            >,
            matrix &
        >::type
        operator+= (matrix<T> const & rhs)
            {
                typedef matrix<T> other_matrix;
                typedef fusion::vector<matrix &, other_matrix const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, rhs), detail::plus_assign()
                );
                return *this;
            }

        /** Subtracts the values in @c rhs from @c *this,
            element-by-element. @c *this and @c m must have the same
            dimensions, and every difference <c>(*this)(i, j) - rhs(i, j)</c>
            must be convertible to the type of <c>(*this)(i, j)</c>. */
        template <typename T>
        typename enable_if<
            mpl::and_<
                mpl::equal_to<typename matrix<T>::num_rows_t, num_rows_t>,
                mpl::equal_to<typename matrix<T>::num_columns_t, num_columns_t>
            >,
            matrix &
        >::type
        operator-= (matrix<T> const & rhs)
            {
                typedef matrix<T> other_matrix;
                typedef fusion::vector<matrix &, other_matrix const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, rhs), detail::minus_assign()
                );
                return *this;
            }

        /** Multiplies each element of @c *this by @c val.  Every product
            <c>(*this)(i, j) * val</c> must be convertible to the type of
            <c>(*this)(i, j)</c>. */
        template <typename T>
        matrix & operator*= (T const & val)
            {
                typedef fusion::vector<matrix &, matrix const &, T const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, *this, val), detail::matrix_scalar_mul_assign()
                );
                return *this;
            }

        /** Divides each element of @c *this by @c val.  Every quotient
            <c>(*this)(i, j) / val</c> must be convertible to the type of
            <c>(*this)(i, j)</c>. */
        template <typename T>
        matrix & operator/= (T const & val)
            {
                typedef fusion::vector<matrix &, matrix const &, T const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, *this, val), detail::matrix_scalar_div_assign()
                );
                return *this;
            }

    private:
        template <typename T>
        friend class matrix;

#ifndef BOOST_UNITS_BLAS_DOXYGEN
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

        // If you're seeing an error at one of these assertions, it's because
        // length-0 rows and columns are illegal.
        BOOST_MPL_ASSERT((mpl::less<mpl::size_t<0>, num_rows_t>));
        BOOST_MPL_ASSERT((mpl::less<mpl::size_t<0>, num_columns_t>));

        // If you're seeing an error here, it's because you used some type
        // sequence besides boost::fusion::vector in your Rows template
        // parameter.  Rows is supposed to be a boost::fusion::vector of
        // boost::fusion::vectors *only*.  If you want to use something else for
        // Rows, you must use the metafunction make_matrix<Rows> instead of
        // supplying the Rows parameter directly to this template.
        BOOST_MPL_ASSERT((is_same<Rows, value_types>));

        template <typename T>
        void assign_data (T const & rhs)
            {
                typedef fusion::vector<matrix &, T const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, rhs), detail::assign()
                );
            }
#endif

        value_types data_;
    };

    /** Returns a @c Matrix whose diagonal elements are 1, and whose
        nondiagonal elements are 0. */
    template <typename Matrix>
    typename enable_if<
        mpl::equal_to<
            typename Matrix::num_rows_t,
            typename Matrix::num_columns_t
        >,
        Matrix
    >::type
    make_identity_matrix ()
    {
        Matrix retval;
        typedef fusion::vector<Matrix &> ops;
        iterate<typename Matrix::num_rows_t>(
            ops(retval), detail::identity_assign()
        );
        return retval;
    }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_MATRIX_HPP
