// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_ITERATE_HPP
#define BOOST_UNITS_BLAS_ITERATE_HPP

#include <boost/fusion/sequence/intrinsic/value_at.hpp>
#include <boost/type_traits/remove_reference.hpp>


namespace boost { namespace units_blas {

    template <std::size_t N>
    struct iterate_unrolled
    {
        template <typename Operands, typename F>
        static void call (Operands const & operands, F const & f)
            {
                F::template call<N - 1>(operands);
                F::template call<N - 2>(operands);
                F::template call<N - 3>(operands);
                F::template call<N - 4>(operands);
                iterate_unrolled<N - 4>::call(operands, f);
            }
    };

    template <>
    struct iterate_unrolled<3>
    {
        template <typename Operands, typename F>
        static void call (Operands const & operands, F const & f)
            {
                F::template call<2>(operands);
                F::template call<1>(operands);
                F::template call<0>(operands);
            }
    };

    template <>
    struct iterate_unrolled<2>
    {
        template <typename Operands, typename F>
        static void call (Operands const & operands, F const & f)
            {
                F::template call<1>(operands);
                F::template call<0>(operands);
            }
    };

    template <>
    struct iterate_unrolled<1>
    {
        template <typename Operands, typename F>
        static void call (Operands const & operands, F const & f)
            { F::template call<0>(operands); }
    };

    template <>
    struct iterate_unrolled<0>
    {
        template <typename Operands, typename F>
        static void call (Operands const &, F const &) {}
    };

    template <typename N, typename Operands, typename F>
    void iterate (Operands const & operands, F const & f)
    { iterate_unrolled<N::value>::call(operands, f); }

    namespace detail {

        template <typename Operands, std::size_t N>
        struct operand_n
        {
            typedef typename remove_reference<
                typename fusion::result_of::value_at_c<Operands, N>::type
            >::type type;
        };

    } // namespace detail

} } // namespace boost::units_blas

#endif // BOOST_UNITS_BLAS_ITERATE_HPP
