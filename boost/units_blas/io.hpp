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

        template <typename Matrix>
        struct print
        {
            template <std::size_t I>
            void call ()
            {
                if (I && I % Matrix::num_columns == 0)
                    os_ << '\n';
                os_ << tuple_access::get<I>(m_) << ' ';
            }

            Matrix m_;
            std::ostream & os_;
        };

    } // namespace detail

    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    std::ostream &
    matrix_t<Tuple, Rows, Columns>::print (std::ostream & os) const
    {
        hana::fold.left(
            data_,
            hana::size_t<0>,
            [&](auto i, auto x) {
                if (i && i % num_columns == 0)
                    os << '\n';
                os << hana::at(data_, i);
                return hana::succ(i);
            }
        );
#if 0
        detail::iterate_simple<num_elements>(
            detail::print<self_type>{*this, os}
        );
#endif
        return os;
    }

    /** Writes matrix<> @c m to stream @c os, using @c matrix<>::print(). */
    template <typename Tuple, std::size_t Rows, std::size_t Columns>
    std::ostream &
    operator<< (std::ostream & os, matrix_t<Tuple, Rows, Columns> m)
    {
        m.print(os);
        return os;
    }

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_IO_HPP
