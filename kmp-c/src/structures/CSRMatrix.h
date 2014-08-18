#pragma once
#include "util/aligned.h"
#include "solvers/VectorOperator.h"
#include "mkl.h"

struct CSRMatrix : public VectorOperator
{
    using VectorOperator::apply;
    virtual void apply(double *x, double *y)
    {
        char n = 'n';
        mkl_cspblas_dcsrgemv(&n, &dim, data.data(), cols.data(), rows.data(), x, y);
    }

    virtual int getDimension() const
    {
        return dim;
    }

    aligned_vector<double> data;
    aligned_vector<int> cols;
    aligned_vector<int> rows;

    int dim;
    int nnz;

};


