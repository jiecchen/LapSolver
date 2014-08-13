/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Graph class for easy use with BLAS
 */
#pragma once
#include <util/aligned.h>
#include "edgelist.h"

typedef struct CSRMatrix
{
    aligned_vector<double> data;
    aligned_vector<int> cols;
    aligned_vector<int> rows;

    int dim;
    int nnz;
} CSRMatrix;

class Graph
{
public:
    explicit Graph(EdgeList &edges);
    const double * getWeights(int v) const;
    const int * getNeighbors(int v) const;
    int getDegree(int v) const;

    void debugPrint();

    const int nv;

private:
    CSRMatrix adj;
    aligned_vector<int> degrees;
};
