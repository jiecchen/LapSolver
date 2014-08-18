#include <mkl.h>
#include <cstdio>
#include "ConjugateGradientSolver.h"

ConjugateGradientSolver::ConjugateGradientSolver(VectorOperator *a, VectorOperator *m, int maxIters, double tolerance)
    : a(a), m(m), maxIters(maxIters), tolerance(tolerance), dim(a->getDimension())
{
    if (m && dim != m->getDimension())
        throw "preconditioner-matrix dimension mismatch";

    tmp = aligned_vector<double>(4 * dim);
}

#define dcg_throw(f) throw fprintf(stderr, "Error! " f " failed with code %d\n", rci), rci;

void ConjugateGradientSolver::apply(double *b, double *x)
{
    // Necessary variables for RCI solver
    int rci;
    double norm = tolerance + 1;

    int one_i = 1;
    double one_d = -1.e0;

    aligned_vector<double> scratch(dim);

    // Initialize the solution
    x[0:dim] = 0.0;

    dcg_init(&dim, x, b, &rci, ipar, dpar, tmp.data());
    if (rci != 0)
        dcg_throw("dcg_init");

    ipar[4] = maxIters;
    ipar[10] = (m != nullptr); // 0 => no preconditioner; else => precondition

    dcg_check(&dim, x, b, &rci, ipar, dpar, tmp.data());
    if (rci != 0)
        dcg_throw("dcg_check");

    do
    {
        dcg (&dim, x, b, &rci, ipar, dpar, tmp.data());
        switch (rci)
        {
        case 0: // PCG is done; just quit
            norm = 0.0;
            break;
        case 1: // Apply the operator to temp(0:n), put the result in tmp(n+1:n)
            a->apply(&tmp[0], &tmp[dim]);
            break;
        case 2: // Compute the norm - decide whether to stop early
            a->apply(x, scratch.data());
            daxpy (&dim, &one_d, b, &one_i, scratch.data(), &one_i);
            norm = dnrm2(&dim, scratch.data(), &one_i);
            break;
        case 3: // Apply the preconditioner, if any. Will never reach if no precon.
            m->apply(&tmp[2 * dim], &tmp[3 * dim]);
            break;
        default:
            dcg_throw("dcg");
        }
    }
    while (norm > tolerance);

    dcg_get (&dim, x, b, &rci, ipar, dpar, tmp.data(), &nIter);
}
