// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_DETAIL_MATRIX_PRODUCT_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_DETAIL_MATRIX_PRODUCT_HPP

#include <boost/mpl/range_c.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/units/quantity.hpp>


namespace boost { namespace units_blas { namespace result_of { namespace detail {

    template <typename T0, typename T1>
    struct value_product
    { typedef BOOST_TYPEOF((T0() * T1())) type; };

    template <typename T0, typename T1>
    struct value_quotient
    { typedef BOOST_TYPEOF((T0(1) / T1(1))) type; };

    template <typename Unit0, typename ValueType0, typename Unit1, typename ValueType1>
    struct value_quotient<
        units::quantity<Unit0, ValueType0>,
        units::quantity<Unit1, ValueType1>
    >
    {
        typedef BOOST_TYPEOF((
            units::quantity<Unit0, ValueType0>::from_value(static_cast<ValueType0>(1)) /
            units::quantity<Unit1, ValueType1>::from_value(static_cast<ValueType1>(1))
        )) type;
    };

    template <typename Unit0, typename ValueType0, typename T1>
    struct value_quotient<units::quantity<Unit0, ValueType0>, T1>
    {
        typedef BOOST_TYPEOF((
            units::quantity<Unit0, ValueType0>::from_value(static_cast<ValueType0>(1)) /
            static_cast<ValueType0>(1)
        )) type;
    };

    template <typename T0, typename Unit1, typename ValueType1>
    struct value_quotient<T0, units::quantity<Unit1, ValueType1> >
    {
        typedef BOOST_TYPEOF((
            static_cast<ValueType1>(1) /
            units::quantity<Unit1, ValueType1>::from_value(static_cast<ValueType1>(1))
        )) type;
    };

    template <typename T0, typename T1>
    struct value_sum
    { typedef BOOST_TYPEOF((T0() + T1())) type; };

    template <typename T0, typename T1>
    struct value_difference
    { typedef BOOST_TYPEOF((T0() - T1())) type; };

    template <typename T, bool fundamental = is_arithmetic<T>::value>
    struct value_inverse;

    template <typename T>
    struct value_inverse<T, true>
    { typedef T type; };

    template <typename T>
    struct value_inverse<T, false>
    { typedef BOOST_TYPEOF((1 / T(1))) type; };

    template <typename Unit, typename ValueType>
    struct value_inverse<units::quantity<Unit, ValueType>, false>
    {
        typedef BOOST_TYPEOF((
            static_cast<ValueType>(1) /
            units::quantity<Unit, ValueType>::from_value(static_cast<ValueType>(1))
        )) type;
    };

    template <typename SequenceL, typename SequenceR, typename N>
    struct value_sum_n :
        mpl::lambda<
            value_sum<
                typename fusion::result_of::value_at<SequenceL, N>::type,
                typename fusion::result_of::value_at<SequenceR, N>::type
            >
        >::type
    {};

    template <typename SequenceL, typename SequenceR, typename N>
    struct value_difference_n :
        mpl::lambda<
            value_difference<
                typename fusion::result_of::value_at<SequenceL, N>::type,
                typename fusion::result_of::value_at<SequenceR, N>::type
            >
        >::type
    {};

    template <typename SequenceL, typename SequenceR, typename N>
    struct value_product_n :
        mpl::lambda<
            value_product<
                typename fusion::result_of::value_at<SequenceL, N>::type,
                typename fusion::result_of::value_at<SequenceR, N>::type
            >
        >::type
    {};

    template <typename SequenceL, typename SequenceR, typename N>
    struct value_quotient_n :
        mpl::lambda<
            value_quotient<
                typename fusion::result_of::value_at<SequenceL, N>::type,
                typename fusion::result_of::value_at<SequenceR, N>::type
            >
        >::type
    {};

    template <typename SequenceL, typename SequenceR>
    struct element_sums
    {
        typedef typename mpl::transform_view<
            mpl::range_c<std::size_t, 0, fusion::result_of::size<SequenceL>::value>,
            value_sum_n<SequenceL, SequenceR, mpl::_1>
        >::type type;
    };

    template <typename SequenceL, typename SequenceR>
    struct element_differences
    {
        typedef typename mpl::transform_view<
            mpl::range_c<std::size_t, 0, fusion::result_of::size<SequenceL>::value>,
            value_difference_n<SequenceL, SequenceR, mpl::_1>
        >::type type;
    };

    template <typename SequenceL, typename SequenceR>
    struct element_products
    {
        typedef typename mpl::transform_view<
            mpl::range_c<std::size_t, 0, fusion::result_of::size<SequenceL>::value>,
            value_product_n<SequenceL, SequenceR, mpl::_1>
        >::type type;
    };

    template <typename SequenceL, typename SequenceR>
    struct element_quotients
    {
        typedef typename mpl::transform_view<
            mpl::range_c<std::size_t, 0, fusion::result_of::size<SequenceL>::value>,
            value_quotient_n<SequenceL, SequenceR, mpl::_1>
        >::type type;
    };

    template <typename SequenceL, typename SequenceR>
    struct dot_product
    {
        typedef typename element_products<SequenceL, SequenceR>::type product_types;

        typedef typename mpl::advance_c<
            typename mpl::begin<product_types>::type,
            1
        >::type begin;
        typedef typename mpl::end<product_types>::type end;
        typedef typename mpl::fold<
            mpl::iterator_range<begin, end>,
            typename mpl::front<product_types>::type,
            typename mpl::lambda<value_sum<mpl::_1, mpl::_2> >::type
        >::type type;
    };

    template <typename ValueTypesL, typename TransposedValueTypesR, typename Row, typename Column>
    struct matrix_product_element :
        mpl::lambda<
            dot_product<
                typename fusion::result_of::value_at<ValueTypesL, Row>::type,
                typename fusion::result_of::value_at<TransposedValueTypesR, Column>::type
            >
        >::type
    {};

    template <typename ValueTypesL, typename TransposedValueTypesR, typename NumColumns>
    struct matrix_product_row
    {
        template <typename Row>
        struct apply
        {
            typedef typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, NumColumns::value>,
                    matrix_product_element<ValueTypesL, TransposedValueTypesR, Row, mpl::_1>
                >::type
            >::type type;
        };
    };

    template <typename Matrix, typename T, typename Row, typename Column>
    struct scalar_product_element :
        mpl::lambda<
            value_product<
                typename result_of::value_at<Matrix, Row, Column>::type,
                T
            >
        >::type
    {};

    template <typename Matrix, typename T, typename NumColumns>
    struct scalar_product_row
    {
        template <typename Row>
        struct apply
        {
            typedef typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, NumColumns::value>,
                    scalar_product_element<Matrix, T, Row, mpl::_1>
                >::type
            >::type type;
        };
    };

    template <typename Matrix, typename T, typename Row, typename Column>
    struct scalar_quotient_element :
        mpl::lambda<
            value_quotient<
                typename result_of::value_at<Matrix, Row, Column>::type,
                T
            >
        >::type
    {};

    template <typename Matrix, typename T, typename NumColumns>
    struct scalar_quotient_row
    {
        template <typename Row>
        struct apply
        {
            typedef typename fusion::result_of::as_vector<
                typename mpl::transform_view<
                    mpl::range_c<std::size_t, 0, NumColumns::value>,
                    scalar_quotient_element<Matrix, T, Row, mpl::_1>
                >::type
            >::type type;
        };
    };

    template <typename ValueTypesL, typename ValueTypesR>
    struct matrix_element_sum_row
    {
        template <typename Row>
        struct apply
        {
            typedef typename fusion::result_of::as_vector<
                typename element_sums<
                    typename fusion::result_of::value_at<ValueTypesL, Row>::type,
                    typename fusion::result_of::value_at<ValueTypesR, Row>::type
                >::type
            >::type type;
        };
    };

    template <typename ValueTypesL, typename ValueTypesR>
    struct matrix_element_difference_row
    {
        template <typename Row>
        struct apply
        {
            typedef typename fusion::result_of::as_vector<
                typename element_differences<
                    typename fusion::result_of::value_at<ValueTypesL, Row>::type,
                    typename fusion::result_of::value_at<ValueTypesR, Row>::type
                >::type
            >::type type;
        };
    };

    template <typename ValueTypesL, typename ValueTypesR>
    struct matrix_element_product_row
    {
        template <typename Row>
        struct apply
        {
            typedef typename fusion::result_of::as_vector<
                typename element_products<
                    typename fusion::result_of::value_at<ValueTypesL, Row>::type,
                    typename fusion::result_of::value_at<ValueTypesR, Row>::type
                >::type
            >::type type;
        };
    };

    template <typename ValueTypesL, typename ValueTypesR>
    struct matrix_element_quotient_row
    {
        template <typename Row>
        struct apply
        {
            typedef typename fusion::result_of::as_vector<
                typename element_quotients<
                    typename fusion::result_of::value_at<ValueTypesL, Row>::type,
                    typename fusion::result_of::value_at<ValueTypesR, Row>::type
                >::type
            >::type type;
        };
    };

} } } } // namespace boost::units_blas::result_of::detail

#endif // BOOST_UNITS_BLAS_RESULT_OF_DETAIL_MATRIX_PRODUCT_HPP
