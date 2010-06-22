// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/units_blas/dont_care.hpp>


namespace bub = boost::units_blas;

int test_main (int, char *[])
{
    bub::_ _;
    _ = 1;

    return 0;
}
