// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DETAIL_LU_HPP
#define BOOST_UNITS_BLAS_DETAIL_LU_HPP

#include <boost/units_blas/exception.hpp>
#include <boost/units_blas/detail/abs.hpp>

#include <boost/throw_exception.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/type_traits/is_fundamental.hpp>
#include <boost/units/is_dimensionless.hpp>


namespace boost { namespace units_blas { namespace detail {

    template <typename T>
    struct is_dimensionless :
        mpl::and_<
            is_fundamental<T>,
            is_arithmetic<T>
        >
    {};

    template <typename Unit, typename ValueType>
    struct is_dimensionless<units::quantity<Unit, ValueType> > :
        units::is_dimensionless<units::quantity<Unit, ValueType> >
    {};

    template <typename T>
    struct unit_of : mpl::void_ {};

    template <typename Unit, typename ValueType>
    struct unit_of<units::quantity<Unit, ValueType> >
    { typedef Unit type; };

    template <typename T>
    struct dimension_of : mpl::void_ {};

    template <typename Unit, typename ValueType>
    struct dimension_of<units::quantity<Unit, ValueType> >
    { typedef typename Unit::dimension_type type; };

    template <typename Matrix, typename T>
    struct nth_type :
        mpl::lambda<
            result_of::value_at<
                Matrix,
                mpl::divides<T, typename Matrix::num_columns_t>,
                mpl::modulus<T, typename Matrix::num_columns_t>
            >
        >::type
    {};

    template <typename Pair, typename ValueType>
    struct update_set_pair :
        mpl::lambda<
            mpl::pair<
                typename mpl::insert<
                    typename Pair::first,
                    typename dimension_of<ValueType>::type
                >::type,
                typename mpl::insert<
                    typename Pair::second,
                    typename unit_of<ValueType>::type
                >::type
            >
        >::type
    {};

    template <typename Matrix>
    struct is_lu_decomposable
    {
        typedef typename mpl::fold<
            mpl::range_c<
                std::size_t,
                0,
                mpl::times<
                    typename Matrix::num_rows_t,
                    typename Matrix::num_columns_t
                >::type::value
            >,
            mpl::pair<mpl::set<>, mpl::set<> >,
            mpl::eval_if<
                is_dimensionless<nth_type<Matrix, mpl::_2> >,
                mpl::identity<mpl::_1>,
                update_set_pair<mpl::_1, nth_type<Matrix, mpl::_2> >
            >
        >::type units_and_dimensions;

        typedef typename mpl::equal_to<
            mpl::size<typename units_and_dimensions::first>,
            mpl::size<typename units_and_dimensions::second>
        >::type type;
    };

    template <typename TempMatrix>
    typename TempMatrix::value_type::value_type
    lu_decompose (TempMatrix & m,
                  array<std::size_t, TempMatrix::static_size> & indices)
    {
        BOOST_STATIC_ASSERT((
            static_cast<std::size_t>(TempMatrix::static_size) ==
            static_cast<std::size_t>(TempMatrix::value_type::static_size)
        ));

        typedef typename TempMatrix::value_type::value_type value_type;

        value_type const epsilon(1.0e-20);
        std::size_t const N = TempMatrix::static_size;

        value_type retval(1.0);
        array<value_type, N> scale_factors;
        for (std::size_t i = 0; i < N; ++i) {
            value_type max(0.0);
            for (std::size_t j = 0; j < N; ++j) {
                value_type tmp;
                if (max < (tmp = abs_(m[i][j])))
                    max = tmp;
                if (!max)
                    throw_exception(singular_matrix());
                scale_factors[i] = 1.0 / max;
            }
        }

        for (std::size_t j = 0; j < N; ++j) {
            for (std::size_t i = 0; i < j; ++i) {
                value_type sum = m[i][j];
                for (std::size_t k = 0; k < i; ++k) {
                    sum -= m[i][k] * m[k][j];
                    m[i][j] = sum;
                }
            }

            value_type max(0.0);
            std::size_t max_index;
            for (std::size_t i = j; i < N; ++i) {
                value_type sum = m[i][j];
                for (std::size_t k = 0; k < j; ++k) {
                    sum -= m[i][k] * m[k][j];
                }
                m[i][j] = sum;
                value_type tmp;
                if (max <= (tmp = scale_factors[i] * abs_(sum))) {
                    max = tmp;
                    max_index = i;
                }
            }

            if (j != max_index) {
                for (std::size_t k = 0; k < N; ++k) {
                    value_type tmp;
                    tmp = m[j][k];
                    m[j][k] = m[max_index][k];
                    m[max_index][k] = tmp;
                }
                retval = -retval;
                scale_factors[max_index] = scale_factors[j];
            }

            indices[j] = max_index;
            if (!m[j][j])
                m[j][j] = epsilon;

            if (j != N - 1) {
                value_type tmp = value_type(1.0) / m[j][j];
                for (std::size_t i = j + 1; i < N; ++i) {
                    m[i][j] *= tmp;
                }
            }
        }

        return retval;
    }

    template <typename TempMatrix>
    void lu_substitute (TempMatrix const & m,
                        array<std::size_t, TempMatrix::static_size> const & indices,
                        array<
                            typename TempMatrix::value_type::value_type,
                            TempMatrix::static_size
                        > & columns)
    {
        BOOST_STATIC_ASSERT((
            static_cast<std::size_t>(TempMatrix::static_size) ==
            static_cast<std::size_t>(TempMatrix::value_type::static_size)
        ));
        typedef typename TempMatrix::value_type::value_type value_type;

        std::size_t const N = TempMatrix::static_size;
        std::size_t ii = 0;
        for (std::size_t i = 0; i < N; ++i) {
            std::size_t index = indices[i];
            value_type sum = columns[index];
            columns[index] = columns[i];
            if (ii) {
                for (std::size_t j = ii - 1; j < i; ++j) {
                    sum -= m[i][j] * columns[j];
                }
            } else if (sum) {
                ii = i + 1;
            }
            columns[i] = sum;
        }

        for (std::size_t i = N - 1; i < N; --i) {
            value_type sum = columns[i];
            for (std::size_t j = i + 1; j < N; ++j) {
                sum -= m[i][j] * columns[j];
            }
            columns[i] = sum / m[i][i];
        }
    }

} } } // namespace boost::units_blas::detail

#endif // BOOST_UNITS_BLAS_DETAIL_LU_HPP
