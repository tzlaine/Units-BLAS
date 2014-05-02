#include <libs/units_blas/test/operations_tests.hpp>

#include <boost/units_blas.hpp>

using namespace boost;

using length_squared = length_sq;

int main()
{
    {
//[default_matrix_decl
    units_blas::matrix<
        std::tuple<length>,
        std::tuple<velocity>,
        std::tuple<length>,
        std::tuple<velocity>
    > my_matrix;
//]
    }

    {
//[make_matrix_decl
    units_blas::matrix<
        std::tuple<length, length>,
        std::tuple<length, length>
    > my_matrix;
//]
    }

    {
//[uniform_matrix_decl
    units_blas::uniform_matrix<time_, 2, 3> my_matrix;
//]
    }

    {
//[equiv_uniform_matrix_decl
    units_blas::matrix<
        std::tuple<time_, time_, time_>,
        std::tuple<time_, time_, time_>
    > my_matrix;
//]
    }

    {
//[vector_decl
    units_blas::vector<length, time_> my_vector;
//]
    }

    {
//[equiv_vector_decl
    units_blas::matrix<
        std::tuple<length>,
        std::tuple<time_>
    > my_vector;
//]
    }

    {
//[transpose_vector_decl
    units_blas::transpose_vector<length, time_> my_transpose_vector;
//]
    }

    {
//[equiv_transpose_vector_decl
    units_blas::matrix<
        std::tuple<length, time_>
    > my_transpose_vector;
//]
    }

    {
//[uniform_vector_decls
    units_blas::uniform_vector<length, 2> my_vector;
    units_blas::uniform_transpose_vector<time_, 2> my_transpose_vector;
//]
    }

    {
//[equiv_uniform_vector_decls
    units_blas::matrix<
        std::tuple<length>,
        std::tuple<length>
    > my_vector;
    
    units_blas::matrix<
        std::tuple<time_, time_>
    > my_transpose_vector;
//]
    }

    {
//[cross_product_space_example
    units_blas::uniform_vector<length, 3> vec;
    units_blas::uniform_vector<length_squared, 3> cross_product_result = vec ^ vec;
//]
    }

    {
//[slicing_example
    using matrix = units_blas::matrix<
        std::tuple<double, double>,
        std::tuple<double, double>
    >;

    matrix m;
    m.at<0, 0>() = 1.0;
    m.at<0, 1>() = 2.0;
    m.at<1, 0>() = 3.0;
    m.at<1, 1>() = 4.0;

    // Note that Rows reorders rows 0 and 1.
    using rows = std::index_sequence<1, 0>;
    using columns = std::index_sequence<1>;

    auto sliced_m =
        units_blas::slice<rows, columns>(m);
    assert((sliced_m.at<0, 0>() == 4.0));
    assert((sliced_m.at<1, 0>() == 2.0));
//]
    }

    {
//[impossible_to_generate_identity_matrix_example
    // No identity matrix exists which will preserve these types.
    units_blas::matrix<
        std::tuple<length, length>,
        std::tuple<length, time_>
    > who_am_i;
//]
    }

    return 0;
}
