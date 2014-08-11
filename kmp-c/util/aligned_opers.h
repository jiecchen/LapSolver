#pragma once
#include "aligned_types.h"
#include "mkl.h"

/* Operator overloads for aligned_vector type */
inline aligned_vector<float> operator+(const aligned_vector<float> &x, const aligned_vector<float> &y)
{
    aligned_vector<float> z(x.size());
    vsAdd(x.size(), x.data(), y.data(), z.data());
    return z;
}

inline aligned_vector<double> operator+(const aligned_vector<double> &x, const aligned_vector<double> &y)
{
    aligned_vector<double> z(x.size());
    vdAdd(x.size(), x.data(), y.data(), z.data());
    return z;
}

inline aligned_vector<float> operator-(const aligned_vector<float> &x, const aligned_vector<float> &y)
{
    aligned_vector<float> z(x.size());
    vsSub(x.size(), x.data(), y.data(), z.data());
    return z;
}

inline aligned_vector<double> operator-(const aligned_vector<double> &x, const aligned_vector<double> &y)
{
    aligned_vector<double> z(x.size());
    vdSub(x.size(), x.data(), y.data(), z.data());
    return z;
}

inline aligned_vector<float> operator*(const aligned_vector<float> &x, const aligned_vector<float> &y)
{
    aligned_vector<float> z(x.size());
    vsMul(x.size(), x.data(), y.data(), z.data());
    return z;
}

inline aligned_vector<double> operator*(const aligned_vector<double> &x, const aligned_vector<double> &y)
{
    aligned_vector<double> z(x.size());
    vdMul(x.size(), x.data(), y.data(), z.data());
    return z;
}

