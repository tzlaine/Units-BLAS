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
#include <boost/mpl/times.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/fusion/algorithm/query/find_if.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/mpl.hpp>
#include <boost/fusion/sequence/comparison/equal_to.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/fusion/sequence/intrinsic/value_at.hpp>

#include <iosfwd>


namespace boost { namespace units_blas {

    template <BOOST_PP_ENUM_PARAMS(FUSION_MAX_VECTOR_SIZE, typename T)>
    struct row :
        fusion::vector<BOOST_PP_ENUM_PARAMS(FUSION_MAX_VECTOR_SIZE, T)>
    {};

    template <BOOST_PP_ENUM_PARAMS(FUSION_MAX_VECTOR_SIZE, typename T)>
    struct all_rows :
        fusion::vector<BOOST_PP_ENUM_PARAMS(FUSION_MAX_VECTOR_SIZE, T)>
    {};

    namespace detail {
        template <typename Rows>
        struct deep_as_vector
        {
            typedef typename fusion::result_of::as_vector<
                typename mpl::transform<
                    Rows,
                    fusion::result_of::as_vector<mpl::_1>
                >::type
            >::type type;
        };
    }

    template <typename Rows>
    class matrix
    {
    public:
        typedef std::size_t size_type;
        typedef typename detail::deep_as_vector<Rows>::type value_types;
        typedef mpl::true_ is_matrix;

        typedef size_t_<
            fusion::result_of::size<value_types>::type::value
        > num_rows_t;
        typedef size_t_<
            fusion::result_of::size<
                typename fusion::result_of::value_at_c<
                    value_types,
                    0
                >::type
            >::type::value
        > num_columns_t;

        matrix () :
            data_ ()
            {}

        matrix (matrix const & m) :
            data_ (m.data_)
            {}

        template <typename T>
        matrix (matrix<T> const & m) :
            data_ ()
            {
                typedef matrix<T> other_matrix;
                BOOST_MPL_ASSERT((mpl::equal_to<typename other_matrix::num_rows_t, num_rows_t>));
                BOOST_MPL_ASSERT((mpl::equal_to<typename other_matrix::num_columns_t, num_columns_t>));
                assign_data(m);
            }

        matrix & operator= (matrix const & rhs)
            {
                data_ = rhs.data_;
                return *this;
            }

        template <typename T>
        matrix & operator= (matrix<T> const & rhs)
            {
                typedef matrix<T> other_matrix;
                BOOST_MPL_ASSERT((mpl::equal_to<typename other_matrix::num_rows_t, num_rows_t>));
                BOOST_MPL_ASSERT((mpl::equal_to<typename other_matrix::num_columns_t, num_columns_t>));
                assign_data(rhs);
                return *this;
            }

        size_type size () const
            { return mpl::times<num_rows_t, num_columns_t>::value; }

        size_type rows () const
            { return num_rows_t::value; }

        size_type columns () const
            { return num_columns_t::value; }

        template <size_type I, size_type J>
        typename result_of::at_c<matrix const, I, J>::type
        at () const
            { return fusion::at_c<J>(fusion::at_c<I>(data_)); }

        std::ostream & print (std::ostream & os) const;


        template <size_type I, size_type J>
        typename result_of::at_c<matrix, I, J>::type
        at ()
            { return fusion::at_c<J>(fusion::at_c<I>(data_)); }

        template <typename T>
        matrix & operator+= (matrix<T> const & rhs)
            {
                typedef matrix<T> other_matrix;
                BOOST_MPL_ASSERT((mpl::equal_to<typename other_matrix::num_rows_t, num_rows_t>));
                BOOST_MPL_ASSERT((mpl::equal_to<typename other_matrix::num_columns_t, num_columns_t>));
                typedef fusion::vector<matrix &, other_matrix const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, rhs), detail::plus_assign()
                );
                return *this;
            }

        template <typename T>
        matrix & operator-= (matrix<T> const & rhs)
            {
                typedef matrix<T> other_matrix;
                BOOST_MPL_ASSERT((mpl::equal_to<typename other_matrix::num_rows_t, num_rows_t>));
                BOOST_MPL_ASSERT((mpl::equal_to<typename other_matrix::num_columns_t, num_columns_t>));
                typedef fusion::vector<matrix &, other_matrix const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, rhs), detail::minus_assign()
                );
                return *this;
            }

        template <typename T>
        matrix & operator*= (T const & rhs)
            {
                typedef fusion::vector<matrix &, matrix const &, T const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, *this, rhs), detail::matrix_scalar_mul_assign()
                );
                return *this;
            }

        template <typename T>
        matrix & operator/= (T const & rhs)
            {
                typedef fusion::vector<matrix &, matrix const &, T const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, *this, rhs), detail::matrix_scalar_div_assign()
                );
                return *this;
            }

    private:
        template <typename T>
        friend class matrix;

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

        // If you're seeing an error at this assert, it's because you defined
        // two or more rows that have different numbers of columns.
        BOOST_MPL_ASSERT((column_lengths_equal));

        // If you're seeing an error at one of these asserts, it's because
        // length-0 rows and columns are illegal.
        BOOST_MPL_ASSERT((mpl::less<size_t_<0>, num_rows_t>));
        BOOST_MPL_ASSERT((mpl::less<size_t_<0>, num_columns_t>));

#if BOOST_UNITS_BLAS_ENFORCE_FUSION_VECTOR_MATRIX_TEMPLATE_PARAMETERS
        // If you're seeing an error here, it's because someone defined the
        // above macro to be nonzero, and you used some type sequence besides
        // boost::fusion::vector in your Rows template parameter.
        BOOST_MPL_ASSERT((is_same<Rows, value_types>));
#endif

        template <typename T>
        void assign_data (T const & rhs)
            {
                typedef fusion::vector<matrix &, T const &> ops;
                iterate<mpl::times<num_rows_t, num_columns_t> >(
                    ops(*this, rhs), detail::assign()
                );
            }

        value_types data_;
    };

    template <typename Matrix>
    typename enable_if<
        mpl::equal_to<
            typename Matrix::num_rows_t,
            typename Matrix::num_columns_t
        >,
        Matrix
    >::type
    identity_matrix ()
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
