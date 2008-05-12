#include <libs/units_blas/test/operations_tests.hpp>

#include <boost/units_blas.hpp>
#include <boost/fusion/include/set.hpp>
#include <boost/fusion/include/list.hpp>
#include <boost/mpl/vector_c.hpp>

using namespace boost;

typedef length Length;
typedef time_ Time;
typedef velocity Velocity;
typedef length_sq LengthSquared;

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

    {
//[cross_product_space_example
    units_blas::uniform_vector<Length, 3>::type vec;
    units_blas::uniform_vector<LengthSquared, 3>::type cross_product_result = vec ^ vec;
//]
    }

    {
//[slicing_example
    typedef units_blas::matrix<
        fusion::vector<
            fusion::vector<double, double>,
            fusion::vector<double, double>
        >
    > Matrix;

    Matrix matrix;
    matrix.at<0, 0>() = 1.0;
    matrix.at<0, 1>() = 2.0;
    matrix.at<1, 0>() = 3.0;
    matrix.at<1, 1>() = 4.0;

    // Note that Rows reorders rows 0 and 1.
    typedef mpl::vector_c<std::size_t, 1, 0> Rows;
    typedef mpl::vector_c<std::size_t, 1> Columns;

    typedef units_blas::result_of::slice<Matrix, Rows, Columns>::type SlicedMatrix;
    SlicedMatrix sliced_matrix =
        units_blas::slice<Rows, Columns>(matrix);
    assert((sliced_matrix.at<0, 0>() == 4.0));
    assert((sliced_matrix.at<1, 0>() == 2.0));
//]
    }

    {
//[impossible_to_generate_identity_matrix_example
    // No identity matrix exists which will preserve these types.
    units_blas::matrix<
        fusion::vector<
            fusion::vector<Length, Length>,
            fusion::vector<Length, Time>
        >
    > who_am_i;
//]
    }

    return 0;
}
