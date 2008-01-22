#include <boost/units_blas.hpp>

#include <boost/fusion/include/set.hpp>
#include <boost/fusion/include/list.hpp>

using namespace boost;

struct Length {};
struct Time {};
struct Velocity {};

int main()
{
    {
//[default_matrix_decl
    units_blas::matrix<
        fusion::vector<
            fusion::vector<Length>,
            fusion::vector<Velocity>,
            fusion::vector<Length>,
            fusion::vector<Velocity>
        >
    > my_matrix;
//]
    }

    {
//[make_matrix_decl
    units_blas::make_matrix<
        mpl::set<
            mpl::vector<Length, Length>,
            fusion::list<Length, Length>
        >
    >::type my_matrix;
//]
    }

    {
//[uniform_matrix_decl
    units_blas::uniform_matrix<Time, 2, 3>::type my_matrix;
//]
    }

    {
//[equiv_uniform_matrix_decl
    units_blas::matrix<
        fusion::vector<
            fusion::vector<Time, Time, Time>,
            fusion::vector<Time, Time, Time>
        >
    > my_matrix;
//]
    }

    {
//[vector_decl
    units_blas::vector<fusion::vector<Length, Time> >::type my_vector;
//]
    }

    {
//[equiv_vector_decl
    units_blas::matrix<
        fusion::vector<
            fusion::vector<Length>,
            fusion::vector<Time>
        >
    > my_vector;
//]
    }

    {
//[transpose_vector_decl
    units_blas::transpose_vector<fusion::vector<Length, Time> >::type my_transpose_vector;
//]
    }

    {
//[equiv_transpose_vector_decl
    units_blas::matrix<
        fusion::vector<
            fusion::vector<Length, Time>
        >
    > my_transpose_vector;
//]
    }

    {
//[uniform_vector_decls
    units_blas::uniform_vector<Length, 2>::type my_vector;
    units_blas::uniform_transpose_vector<Time, 2>::type my_transpose_vector;
//]
    }

    {
//[equiv_uniform_vector_decls
    units_blas::matrix<
        fusion::vector<
            fusion::vector<Length>,
            fusion::vector<Length>
        >
    > my_vector;
    
    units_blas::matrix<
        fusion::vector<
            fusion::vector<Time, Time>
        >
    > my_transpose_vector;
//]
    }

    return 0;
}
