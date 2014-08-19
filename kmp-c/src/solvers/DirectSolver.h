#pragma once
#include <memory>
#include "structures/CSRMatrix.h"
#include "VectorOperator.h"
#include "mkl.h"

#define dss_guard(f, ...) \
    err = f(handle, ## __VA_ARGS__); \
    if (err != MKL_DSS_SUCCESS) \
        throw fprintf(stderr, "Error! " #f " failed with code %d\n", err), err;

template <MatrixType M>
class DirectSolver : public VectorOperator
{
public:
    using VectorOperator::apply;
    DirectSolver(CSRMatrix<M> *csr)
        : csr(csr),
          dim(csr->getDimension())
    {
        sym = (M == Real || M == RealPosDef)
              ? MKL_DSS_NON_SYMMETRIC
              : MKL_DSS_SYMMETRIC;
        type = (M == Real || M == RealSym)
               ? MKL_DSS_INDEFINITE
               : MKL_DSS_POSITIVE_DEFINITE;

        opt = MKL_DSS_DEFAULTS + MKL_DSS_ZERO_BASED_INDEXING;
        dss_guard(dss_create, opt);
        opt = MKL_DSS_DEFAULTS;
        dss_guard(dss_define_structure, sym,
                  csr->rows.data(), csr->dim, csr->dim,
                  csr->cols.data(), csr->nnz);
        dss_guard(dss_reorder, opt, NULL);
        dss_guard(dss_factor_real, type, csr->data.data());
    }

    virtual ~DirectSolver()
    {
        dss_guard(dss_delete, opt);
    }

    virtual int getDimension() const
    {
        return dim;
    }

    virtual void apply(double *b, double *x)
    {
        MKL_INT nRhs = 1;
        dss_guard(dss_solve_real, opt, b, nRhs, x);
    }

private:
    int dim;
    CSRMatrix<M> *csr;

    MKL_INT opt;
    MKL_INT sym;
    MKL_INT type;
    _MKL_DSS_HANDLE_t handle;
    _INTEGER_t err;
};
