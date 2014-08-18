#pragma once
#include "VectorOperator.h"

class ConjugateGradientSolver : public VectorOperator
{
public:
	using VectorOperator::apply;
    ConjugateGradientSolver(VectorOperator *a, VectorOperator *m = nullptr,
                            int maxIters = 1000, double tolerance = 1e-8);
    virtual ~ConjugateGradientSolver() {}

    virtual int getDimension() const
    {
        return dim;
    }

    virtual int getIterations() const
    {
        return nIter;
    }

    virtual void apply(double *b, double *x);

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