#include <solvers/ConjugateGradientSolver.h>
#include <util/aligned.h>
#include <util/Benchmark.h>
#include <structures/CSRMatrix.h>
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

class JacobiPreconditioner : public VectorOperator
{
public:
    JacobiPreconditioner(aligned_vector<int> &rows, aligned_vector<int> &cols, aligned_vector<double> &data)
        : n(rows.size()-1), data(data.data()), cols(cols.data()), rows(rows.data())
    {

    }

    void apply(double *x, double *y)
    {
        double one = 1.E0;
        char matdes[] = { 'd', 'l', 'n', 'c' };
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

    aligned_vector<int> ia = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
    aligned_vector<int> ja = { 1, 3, 6, 7, 2, 3, 5, 3, 8, 4, 7, 5, 6, 7, 6, 8, 7, 8 };
    aligned_vector<double> da = { 7, 1, 2, 7, -4, 8, 2, 1, 5, 7, 9, 5, 1, 5, -1, 5, 11, 5 };

    // Change to 0-based indexing
    ia.data()[0:ia.size()]--;
    ja.data()[0:ja.size()]--;

    CSRMatrix<RealSym> a(ia, ja, da);
    JacobiPreconditioner precon(ia, ja, da);
    ConjugateGradientSolver cg(&a, &precon);

    // initialize solution via original operator
    double expected_sol[8] = { 1, 0, 1, 0, 1, 0, 1, 0 };
    aligned_vector<double> rhs = a * expected_sol;

    aligned_vector<double> x(8);
    auto bench = make_benchmark(argc, argv, [&] () {
        x = cg.apply(rhs);
    });

    printf("iters = %d\n", cg.getIterations());
    return 0;
}
