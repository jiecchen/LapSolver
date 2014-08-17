#pragma once
#include "VectorOperator.h"

class ConjugateGradientSolver : public VectorOperator
{
public:
    ConjugateGradientSolver(VectorOperator *a, VectorOperator *m,
                            int maxIters = 1000, double tolerance = 1e-10);

    virtual int getDimension() const
    {
        return dim;
    }

    virtual int getIterations() const
    {
        return nIter;
    }

    virtual void apply(double *x, double *b);

private:
    VectorOperator *a;
    VectorOperator *m;
    int maxIters;
    double tolerance;

    // Space needed on every call
    aligned_vector<double> x;
    aligned_vector<double> tmp;
    int ipar[128];
    double dpar[128];

    int dim;
    int nIter = 0;
};