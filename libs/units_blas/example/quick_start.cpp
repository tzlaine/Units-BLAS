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
        fusion::vector1<
            fusion::vector3<Length, Time, Frequency>
        >
    > MatrixType1;
    
    // a 3x1 matrix
    typedef units_blas::matrix<
        fusion::vector3<
            fusion::vector1<Length>,
            fusion::vector1<Time>,
            fusion::vector1<Frequency>
        >
    > MatrixType2;
    
    // a 2x2 matrix
    typedef units_blas::matrix<
        fusion::vector2<
            fusion::vector2<Length, Time>,
            fusion::vector2<Frequency, Time>
        >
    > MatrixType3;
//]

    return 0;
}
