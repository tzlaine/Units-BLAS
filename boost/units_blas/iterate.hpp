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

    /** This struct has one member, @c call().  @c call() calls several
        iterations of @c F::call(), then recursively calls itself, until all
        iterations have been invoked.  Note that the order of iteration is
        descending, from N - 1 to 0.  \see See also @c
        boost::units_blas::iterate(). */
    template <std::size_t N>
    struct iterate_unrolled
    {
        template <typename Operands, typename F>
        static BOOST_UNITS_BLAS_INLINE
        void call (Operands const & operands, F const & f)
            {
                F::template call<N - 1>(operands);
                F::template call<N - 2>(operands);
                F::template call<N - 3>(operands);
                F::template call<N - 4>(operands);
                iterate_unrolled<N - 4>::call(operands, f);
            }
    };

#ifndef BOOST_UNITS_BLAS_DOXYGEN
    template <>
    struct iterate_unrolled<3>
    {
        template <typename Operands, typename F>
        static BOOST_UNITS_BLAS_INLINE
        void call (Operands const & operands, F const & f)
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
        static BOOST_UNITS_BLAS_INLINE
        void call (Operands const & operands, F const & f)
            {
                F::template call<1>(operands);
                F::template call<0>(operands);
            }
    };

    template <>
    struct iterate_unrolled<1>
    {
        template <typename Operands, typename F>
        static BOOST_UNITS_BLAS_INLINE
        void call (Operands const & operands, F const & f)
            { F::template call<0>(operands); }
    };

    template <>
    struct iterate_unrolled<0>
    {
        template <typename Operands, typename F>
        static BOOST_UNITS_BLAS_INLINE
        void call (Operands const &, F const &)
            {}
    };
#endif

    /** This function can be used as a generic mechanism for unrolling loops
        at compile time.  Given a number of iterations @c N, an arbitrary @c
        fusion::vector of operands @c Operands, and a functor type @c F as
        template parameters, it calls @c iterate_unrolled<>::call().  @c
        iterate_unrolled<>::call() calls several iterations of @c F::call(),
        then recursively calls itself, until all iterations have been invoked.
        Note that the order of iteration is descending, from N - 1 to 0. */
    template <typename N, typename Operands, typename F>
    BOOST_UNITS_BLAS_INLINE void iterate (Operands const & operands, F const & f)
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
