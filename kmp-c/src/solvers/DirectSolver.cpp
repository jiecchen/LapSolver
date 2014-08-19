#include "DirectSolver.h"
#include "mkl.h"

template<MatrixType M>
struct DirectSolver<M>::DSS_impl
{
    MKL_INT opt = MKL_DSS_DEFAULTS;
    MKL_INT sym;
    MKL_INT type;
    _MKL_DSS_HANDLE_t handle;
    _INTEGER_t err = MKL_DSS_SUCCESS;
};

constexpr MKL_INT isSym(MatrixType mt)
{
    return (mt == Real || mt == RealPosDef)
           ? MKL_DSS_NON_SYMMETRIC
           : MKL_DSS_SYMMETRIC;
}

constexpr MKL_INT defType(MatrixType mt)
{
    return (mt == Real)
           ? MKL_DSS_INDEFINITE
           : ((mt == RealPosDef)
              ? MKL_DSS_POSITIVE_DEFINITE
              : ((mt == RealSym)
                 ? MKL_DSS_HERMITIAN_INDEFINITE
                 : MKL_DSS_HERMITIAN_POSITIVE_DEFINITE));
}

#define dss_guard(f, ...) \
    _impl->err = f(_impl->handle, __VA_ARGS__); \
    if (_impl->err != MKL_DSS_SUCCESS) \
        throw fprintf(stderr, "Error! " #f " failed with code %d\n", _impl->err), _impl->err;

template<MatrixType M>
DirectSolver<M>::DirectSolver(CSRMatrix<M> *csr)
    : csr(csr),
      dim(csr->getDimension())
{
    _impl = std::unique_ptr<DSS_impl>();
    _impl->sym = isSym(M);
    _impl->type = defType(M);

    dss_guard(dss_create, _impl->opt);
    dss_guard(dss_define_structure, _impl->sym,
              csr->rows.data(), csr->dim, csr->dim,
              csr->cols.data(), csr->nnz);
    dss_guard(dss_reorder, _impl->opt, NULL);
    dss_guard(dss_factor_real, _impl->type, csr->data.data());
}

template<MatrixType M>
DirectSolver<M>::~DirectSolver()
{
    dss_guard(dss_delete, _impl->opt);
}

template<MatrixType M>
void DirectSolver<M>::apply(double *b, double *x)
{
    dss_guard(dss_solve_real, _impl->opt, b, 1, x);
}
