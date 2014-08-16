#pragma once
#include "util/aligned.h"

struct EdgeList
{
    explicit EdgeList(int ne);
    EdgeList(int ne, int *u, int *v);
    EdgeList(int ne, int *u, int *v, double *w);
    EdgeList(const aligned_vector<int> &u, const aligned_vector<int> &v);
    EdgeList(const aligned_vector<int> &u, const aligned_vector<int> &v, const aligned_vector<double> &w);

    const int ne;
    const int nv;

    aligned_vector<int> u;
    aligned_vector<int> v;
    aligned_vector<double> w;
};
