#pragma once
#include "structures/CSRMatrix.h"
#include "VectorOperator.h"

template <MatrixType M>
class DirectSolver : public VectorOperator
{
public:
    using VectorOperator::apply;
    DirectSolver(CSRMatrix<M> *csr);
    virtual ~DirectSolver() {}

    virtual int getDimension() const
    {
        return dim;
    }

    virtual void apply(double *b, double *x);

private:
    CSRMatrix<M> *csr;
    int iparm[64];
    int dim;
};
