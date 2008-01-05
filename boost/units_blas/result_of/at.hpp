// boost.units_blas
//
// Copyright (C) 2008 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNITS_BLAS_RESULT_OF_AT_HPP
#define BOOST_UNITS_BLAS_RESULT_OF_AT_HPP

#include <boost/units_blas/result_of/value_at.hpp>

#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <boost/type_traits/is_const.hpp>


namespace boost { namespace units_blas { namespace result_of {

    template <typename Matrix, typename I, typename J>
    struct at
    {
        typedef typename value_at<Matrix, I, J>::type value_type;
        typedef typename mpl::eval_if<
            is_const<Matrix>,
            add_reference<typename add_const<value_type>::type>,
            add_reference<value_type>
        >::type type;
    };

    template <typename Matrix, std::size_t I, std::size_t J>
    struct at_c
    {
        typedef typename at<Matrix, size_t_<I>, size_t_<J> >::type type;
    };

} } } // namespace boost::units_blas::result_of

#endif // BOOST_UNITS_BLAS_RESULT_OF_AT_HPP
