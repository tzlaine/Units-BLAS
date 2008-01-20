#include <boost/units_blas.hpp>

using namespace boost;

struct Length {};
struct Time {};
struct Frequency {};

int main()
{
// [quick_start_example_types
    // a 1x3 matrix
    typedef units_blas::matrix<
        fusion::vector<
            fusion::vector<Length, Time, Frequency>
        >
    > MatrixType1;
    
    // a 3x1 matrix
    typedef units_blas::matrix<
        fusion::vector<
            fusion::vector<Length>,
            fusion::vector<Time>,
            fusion::vector<Frequency>
        >
    > MatrixType2;
    
    // a 2x2 matrix
    typedef units_blas::matrix<
        fusion::vector<
            fusion::vector<Length, Time>,
            fusion::vector<Frequency, Time>
        >
    > MatrixType3;
//]

    return 0;
}
