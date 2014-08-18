#pragma once

#include "util/aligned.h"

class VectorOperator
{
public:
    virtual ~VectorOperator() {}

    // Non-allocating routines
    virtual void apply(double *, double *) = 0;
    virtual void apply(aligned_vector<double> &x, aligned_vector<double> &y) {
        apply(x.data(), y.data());
    }

    // Allocating routines
    virtual aligned_vector<double> apply(double *x) {
        aligned_vector<double> y(getDimension());
        apply(x, y.data());
        return y;
    }
    virtual aligned_vector<double> apply(aligned_vector<double> &x) {
        return apply(x.data());
    }

    virtual int getDimension() const = 0;
};

template <typename Vec>
aligned_vector<double> operator*(VectorOperator &lhs, Vec rhs)
{
    return lhs.apply(rhs);
}
