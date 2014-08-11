#pragma once
#include "util/aligned_types.h"

class EdgeList
{
public:
    EdgeList(int ne);
    EdgeList(int ne, int *u, int *v);
    EdgeList(int ne, int *u, int *v, double *w);
    EdgeList(const aligned_vector<int> &u, const aligned_vector<int> &v);
    EdgeList(const aligned_vector<int> &u, const aligned_vector<int> &v, const aligned_vector<double> &w);

    bool isSymmetric();

    const int ne;

private:
    aligned_vector<int> u;
    aligned_vector<int> v;
    aligned_vector<double> w;
};