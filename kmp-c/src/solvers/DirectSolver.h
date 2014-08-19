#pragma once
#include <memory>
#include "structures/CSRMatrix.h"
#include "VectorOperator.h"

template <MatrixType M>
class DirectSolver : public VectorOperator
{
public:
    using VectorOperator::apply;
    DirectSolver(CSRMatrix<M> *csr);
    virtual ~DirectSolver();

    virtual int getDimension() const
    {
        return dim;
    }

    virtual void apply(double *b, double *x);

private:
    int dim;
    CSRMatrix<M> *csr;

    // Opaque pointers, aww yeah
    struct DSS_impl;
    std::unique_ptr<DSS_impl> _impl;
};
