#include <solvers/ConjugateGradientSolver.h>
#include <util/aligned.h>
#include <util/Benchmark.h>
#include "mkl.h"

/*
 * Note: this test is identical to Intel's own cg_jacobi_precon.
 *
 * It is here to verify that our interface does not change the
 * correct running of Intel's PCG implementation.
 *
 * It should run in 8 iterations and have one negative zero at
 * dimension 6 (index 5).
 * 
 * This particular preconditioner does not improve the convergence.
 */

// I should probably refactor this into a class in solers
class TestOperator : public VectorOperator
{
public:
    void apply(double *x, double *y)
    {
        mkl_dcsrsymv (&tr, &n, a, ia, ja, x, y);
    }

    int getDimension() const
    {
        return n;
    }

    MKL_INT ia[9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
    MKL_INT ja[18] = { 1, 3, 6, 7, 2, 3, 5, 3, 8, 4, 7, 5, 6, 7, 6, 8, 7, 8 };
    double a[18] = { 7, 1, 2, 7, -4, 8, 2, 1, 5, 7, 9, 5, 1, 5, -1, 5, 11, 5 };

private:
    char tr = 'u';
    int n = 8;
};

class JacobiPreconditioner : public VectorOperator
{
public:
    JacobiPreconditioner(int n, double *data, int *cols, int *rows)
        : n(n), data(data), cols(cols), rows(rows)
    {

    }

    void apply(double *x, double *y)
    {
        double one = 1.E0;
        char matdes[] = { 'd', 'l', 'n', 'f' };
        mkl_dcsrsv (&matdes[2], &n, &one, matdes, data, cols, rows, rows + 1, x, y);
    }

    int getDimension() const
    {
        return n;
    }

private:
    int n;
    double *data;
    int *cols;
    int *rows;
};

int main(int argc, char **argv)
{
    atexit(mkl_free_buffers);
    TestOperator a;
    JacobiPreconditioner precon(8, a.a, a.ja, a.ia);
    ConjugateGradientSolver cg(&a, &precon);

    // initialize solution via original operator
    double expected_sol[8] = { 1, 0, 1, 0, 1, 0, 1, 0 };
    aligned_vector<double> b(8);
    a.apply(expected_sol, b.data());

    aligned_vector<double> x;
    auto bench = make_benchmark(argc, argv, [&] ()
    {
        x = cg.apply(b);
    });

    for (auto xi : x)
    {
        printf("%.3f\n", xi);
    }

    printf("iters = %d\n", cg.getIterations());
    return 0;
}
