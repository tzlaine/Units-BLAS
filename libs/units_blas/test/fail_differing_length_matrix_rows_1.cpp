// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "failure_tests.hpp"

#include <boost/test/minimal.hpp>

namespace bh = boost::hana;

int test_main (int, char *[])
{
    using failure_type = bub::matrix<
        bh::_tuple<double, double>,
        bh::_tuple<double>
    >;

    failure_type f;

    return 0;
}
