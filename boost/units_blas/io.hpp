// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_IO_HPP
#define BOOST_UNITS_BLAS_IO_HPP

#include <boost/units_blas/matrix_fwd.hpp>

#include <boost/fusion/algorithm/iteration/for_each.hpp>

#include <ostream>


namespace boost { namespace units_blas {

    namespace detail {

        struct print_value
        {
            print_value (std::ostream & os) : os_ (os) {}

            template <typename T>
            void operator() (T const & t) const
                { os_ << t << " "; }

            void operator() (_ const &) const
                { os_ << "-- "; }

            std::ostream & os_;
        };

        struct print_row
        {
            print_row (std::ostream & os) : os_ (os) {}

            template <typename T>
            void operator() (T const & r) const
                {
                    fusion::for_each(r, print_value(os_));
                    os_ << "\n";
                }

            std::ostream & os_;
        };

    } // namespace detail

    template <typename T>
    std::ostream & matrix<T>::print (std::ostream & os) const
    {
        fusion::for_each(data_, detail::print_row(os));
        return os;
    }

    /** Writes matrix<> @c m to stream @c os, using @c matrix<>::print(). */
    template <typename T>
    std::ostream & operator<< (std::ostream & os, matrix<T> const & m)
    {
        m.print(os);
        return os;
    }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_IO_HPP
