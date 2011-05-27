// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "failure_tests.hpp"

#include <boost/mpl/vector.hpp>

#include <boost/test/minimal.hpp>


int test_main (int, char *[])
{
    typedef bub::make_matrix<
        boost::mpl::vector<
            boost::mpl::vector<double>
        >
    >::type failure_type;

    failure_type ft;

    return 0;
}
