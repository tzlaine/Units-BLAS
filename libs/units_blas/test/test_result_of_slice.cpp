// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "result_of_tests.hpp"

#include <boost/units_blas/result_of/slice.hpp>

#include <boost/mpl/vector_c.hpp>

#include <boost/test/minimal.hpp>


typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<int, long, float, double>,
        boost::fusion::vector<long, int, float, double>,
        boost::fusion::vector<int, float, long, double>,
        boost::fusion::vector<int, long, double, float>
    >
> matrix_4x4_fundamentals_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<long, int, double>,
        boost::fusion::vector<int, float, double>
    >
> matrix_4x4_fundamentals_type_slice_1;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<long, float>,
        boost::fusion::vector<float, long>
    >
> matrix_4x4_fundamentals_type_slice_2;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<frequency, length, time_, dimensionless>,
        boost::fusion::vector<length, frequency, time_, dimensionless>,
        boost::fusion::vector<frequency, time_, length, dimensionless>,
        boost::fusion::vector<frequency, length, dimensionless, time_>
    >
> matrix_4x4_units_type;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, frequency, dimensionless>,
        boost::fusion::vector<frequency, time_, dimensionless>
    >
> matrix_4x4_units_type_slice_1;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<length, time_>,
        boost::fusion::vector<time_, length>
    >
> matrix_4x4_units_type_slice_2;

int test_main (int, char *[])
{
    typedef boost::mpl::vector_c<std::size_t, 1, 2> rows_1;
    typedef boost::mpl::vector_c<std::size_t, 0, 1, 3> columns_1;
    typedef boost::mpl::vector_c<std::size_t, 0, 2> rows_2;
    typedef boost::mpl::vector_c<std::size_t, 1, 2> columns_2;

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::slice<
                matrix_4x4_fundamentals_type,
                rows_1,
                columns_1
            >::type,
            matrix_4x4_fundamentals_type_slice_1
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::slice<
                matrix_4x4_fundamentals_type,
                rows_2,
                columns_2
            >::type,
            matrix_4x4_fundamentals_type_slice_2
        >::type
    ));

    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::slice<
                matrix_4x4_units_type,
                rows_1,
                columns_1
            >::type,
            matrix_4x4_units_type_slice_1
        >::type
    ));
    BOOST_MPL_ASSERT((
        boost::is_same<
            bub::result_of::slice<
                matrix_4x4_units_type,
                rows_2,
                columns_2
            >::type,
            matrix_4x4_units_type_slice_2
        >::type
    ));

    return 0;
}
