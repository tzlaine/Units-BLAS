// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/matrix.hpp>


namespace bub = boost::units_blas;
namespace bh = boost::hana;

using matrix_3x1 = bub::matrix<
    bh::_tuple<double>,
    bh::_tuple<double>,
    bh::_tuple<double>
>;

using matrix_1x3 = bub::matrix<
    bh::_tuple<double, double, double>
>;

using matrix_2x2 = bub::matrix<
    bh::_tuple<double, double>,
    bh::_tuple<double, double>
>;

using matrix_1x2 = bub::matrix<
    bh::_tuple<double, double>
>;

using matrix_2x1 = bub::matrix<
    bh::_tuple<double>,
    bh::_tuple<double>
>;
