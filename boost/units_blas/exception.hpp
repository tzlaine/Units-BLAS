// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_EXCEPTION_HPP
#define BOOST_UNITS_BLAS_EXCEPTION_HPP

#include <stdexcept>


namespace boost { namespace units_blas {

    /** Thrown when operations involving matrix inversion encounter a singular
        matrix. */
    class singular_matrix :
        public std::exception
    {
    public:
        virtual const char * what () const throw()
        { return "singluar_matrix"; }
    };

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_EXCEPTION_HPP
