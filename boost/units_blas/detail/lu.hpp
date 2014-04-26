// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DETAIL_LU_HPP
#define BOOST_UNITS_BLAS_DETAIL_LU_HPP

#include <boost/mpl/fold.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/set.hpp>
#include <boost/units/is_dimensionless.hpp>

#include <array>
#include <cmath>
#include <cstdlib>


namespace boost { namespace units_blas { namespace detail {

    template <typename T>
    struct is_dimensionless :
        mpl::true_
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

    template <typename T>
    using dimension_of_t = typename dimension_of<T>::type;

    template <typename Matrix, typename T>
    struct nth_type :
        mpl::lambda<
            std::tuple_element<
                T::value,
                typename Matrix::value_types
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
    struct has_uniform_dimensional_units
    {
        typedef typename mpl::fold<
            mpl::range_c<
                std::size_t,
                0,
                Matrix::num_elements
            >,
            mpl::pair<mpl::set<>, mpl::set<> >,
            mpl::eval_if<
                is_dimensionless<nth_type<Matrix, mpl::_2> >,
                mpl::identity<mpl::_1>,
                update_set_pair<mpl::_1, nth_type<Matrix, mpl::_2> >
            >
        >::type units_and_dimensions;

        static bool const value =
            mpl::size<typename units_and_dimensions::first>::value ==
            mpl::size<typename units_and_dimensions::second>::value;
    };

    template <typename TempMatrix>
    auto lu_decompose (
        TempMatrix & m,
        std::array<std::size_t, std::tuple_size<TempMatrix>::value> & indices
    ) {
        static_assert(
            std::tuple_size<TempMatrix>::value ==
            std::tuple_size<typename TempMatrix::value_type>::value,
            "lu_decompse() is only defined for square matrices"
        );

        using value_type = typename TempMatrix::value_type::value_type;
        using std::abs;

        value_type const epsilon(1.0e-20);
        std::size_t const N = std::tuple_size<TempMatrix>::value;

        std::pair<value_type, bool> retval(1.0, true);
        std::array<value_type, N> scale_factors;
        for (std::size_t i = 0; i < N; ++i) {
            value_type max(0.0);
            for (std::size_t j = 0; j < N; ++j) {
                value_type tmp;
                if (max < (tmp = abs(m[i][j])))
                    max = tmp;
                if (!max) {
                    retval.second = false;
                    return retval;
                }
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
                if (max <= (tmp = scale_factors[i] * abs(sum))) {
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
                retval.first = -retval.first;
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
    void lu_substitute (
        TempMatrix m,
        std::array<std::size_t, std::tuple_size<TempMatrix>::value> indices,
        std::array<
            typename TempMatrix::value_type::value_type,
            std::tuple_size<TempMatrix>::value
        > & columns
    ) {
        static_assert(
            std::tuple_size<TempMatrix>::value ==
            std::tuple_size<typename TempMatrix::value_type>::value,
            "lu_decompse() is only defined for square matrices"
        );
        using value_type = typename TempMatrix::value_type::value_type;

        std::size_t const N = std::tuple_size<TempMatrix>::value;
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
