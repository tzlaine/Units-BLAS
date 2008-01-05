// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>


namespace bub = boost::units_blas;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<double>,
        boost::fusion::vector<double>,
        boost::fusion::vector<double>
    >
> matrix_3x1;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<double, double, double>
    >
> matrix_1x3;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<double, double>,
        boost::fusion::vector<double, double>
    >
> matrix_2x2;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<double, double>
    >
> matrix_1x2;

typedef bub::matrix<
    boost::fusion::vector<
        boost::fusion::vector<double>,
        boost::fusion::vector<double>
    >
> matrix_2x1;
