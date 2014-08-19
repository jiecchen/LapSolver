#include <solvers/DirectSolver.h>
#include <util/aligned.h>
#include <util/Benchmark.h>
#include <structures/CSRMatrix.h>
#include "mkl.h"

int main(int argc, char **argv)
{
    atexit(mkl_free_buffers);

    aligned_vector<int> ia = { 1, 6, 7, 8, 9, 10 };
    aligned_vector<int> ja = { 1, 2, 3, 4, 5, 2, 3, 4, 5 };
    aligned_vector<double> da = { 9, 1.5, 6, .75, 3, 0.5, 12, .625, 16 };

    // Change to 0-based indexing
    ia.data()[0:ia.size()]--;
    ja.data()[0:ja.size()]--;

    CSRMatrix<RealSym> a(ia, ja, da);
    DirectSolver<RealSym> dss(&a);

    // initialize solution via original operator
    aligned_vector<double> rhs = { 1, 2, 3, 4, 5 };

    aligned_vector<double> x(8);
    auto bench = make_benchmark(argc, argv, [&] ()
    {
        x = dss.apply(rhs);
    });

    for (auto xi : x)
        printf("%.3f\n", xi);

    return 0;
}
