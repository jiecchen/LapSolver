#pragma once
#include "util/aligned.h"
#include "solvers/VectorOperator.h"
#include "mkl.h"

enum MatrixType
{
    Real,
    RealPosDef,
    RealSym,
    RealSymPosDef,
};

template <MatrixType M = Real>
struct CSRMatrix : public VectorOperator
{
    using VectorOperator::apply;
    CSRMatrix() {}

    CSRMatrix(const aligned_vector<int> &rows, const aligned_vector<int> &cols, const aligned_vector<double> &data)
        : rows(rows), cols(cols), data(data), dim(rows.size() - 1), nnz(cols.size())
    {

    }

    CSRMatrix(aligned_vector<int> &&rows, aligned_vector<int> &&cols, aligned_vector<double> &&data)
        : rows(rows), cols(cols), data(data), dim(rows.size() - 1), nnz(cols.size())
    {

    }

    void apply(double *x, double *y)
    {
        char normal = 'n';
        mkl_cspblas_dcsrgemv(&normal, &dim, data.data(), rows.data(), cols.data(), x, y);
    }

    int getDimension() const
    {
        return dim;
    }

    MatrixType getType() const
    {
        return M;
    }

    aligned_vector<double> data;
    aligned_vector<int> cols;
    aligned_vector<int> rows;

    int dim;
    int nnz;

};

template <>
void CSRMatrix<RealSym>::apply(double *x, double *y)
{
    char upper = 'u';
    mkl_cspblas_dcsrsymv(&upper, &dim, data.data(), rows.data(), cols.data(), x, y);
}

template <>
void CSRMatrix<RealSymPosDef>::apply(double *x, double *y)
{
    char upper = 'u';
    mkl_cspblas_dcsrsymv(&upper, &dim, data.data(), rows.data(), cols.data(), x, y);
}
