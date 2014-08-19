#include "DirectSolver.h"

template<MatrixType M>
DirectSolver<M>::DirectSolver(CSRMatrix<M> *csr)
    : csr(csr),
      dim(csr->getDimension())
{
    
}

template<MatrixType M>
void DirectSolver<M>::apply(double *b, double *x)
{

}
