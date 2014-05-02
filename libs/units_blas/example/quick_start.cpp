#include <boost/units_blas.hpp>

using namespace boost;

struct length {};
struct time_ {};
struct frequency {};

int main()
{
//[quick_start_example_types
    // a 1x3 matrix
    using matrix_type_1 = units_blas::matrix<
        std::tuple<length, time_, frequency>
    >;
    
    // a 3x1 matrix
    using matrix_type_2 = units_blas::matrix<
        std::tuple<length>,
        std::tuple<time_>,
        std::tuple<frequency>
    >;
    
    // a 2x2 matrix
    using matrix_type_3 = units_blas::matrix<
        std::tuple<length, time_>,
        std::tuple<frequency, time_>
    >;
//]

    return 0;
}
