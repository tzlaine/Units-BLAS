// boost.units_blas
//
// Copyright (C) 2010 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_DONT_CARE_HPP
#define BOOST_UNITS_BLAS_DONT_CARE_HPP


namespace boost { namespace units_blas {

    struct _
    {
        _ ()
            {}

        template <typename T>
        explicit _ (T const &)
            {}
    };

    _ operator* (_, _)
    { return _(); }

    template <typename T>
    _ operator* (T, _)
    { return _(); }

    template <typename T>
    _ operator* (_, T)
    { return _(); }

    _ operator/ (_, _)
    { return _(); }

    template <typename T>
    _ operator/ (T, _)
    { return _(); }

    template <typename T>
    _ operator/ (_, T)
    { return _(); }

    _ operator+ (_, _)
    { return _(); }

    template <typename T>
    const T & operator+ (const T & lhs, _)
    { return lhs; }

    template <typename T>
    const T & operator+ (_, const T & rhs)
    { return rhs; }

    _ operator- (_, _)
    { return _(); }

    template <typename T>
    const T & operator- (const T & lhs, _)
    { return lhs; }

    template <typename T>
    const T & operator- (_, const T & rhs)
    { return rhs; }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_DONT_CARE_HPP
