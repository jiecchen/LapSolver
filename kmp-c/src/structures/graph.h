/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Graph class for easy use with BLAS
 */
#pragma once
#include <utility>
#include "util/aligned.h"
#include "EdgeList.h"

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
    Graph();
    Graph(const Graph &rval);
    Graph(Graph &&rval);

    Graph &operator=(const Graph &g)
    {
        if (this != &g)
        {
            nv = g.nv;
            adj = g.adj;
            degrees = g.degrees;
        }
        return *this;
    }

    Graph &operator=(Graph && g)
    {
        nv = std::move(g.nv);
        adj = std::move(g.adj);
        degrees = std::move(g.degrees);
        return *this;
    }

    explicit Graph(EdgeList &&edges);

    const double *getWeights(int v) const;
    const int *getNeighbors(int v) const;
    int getDegree(int v) const;

    void debugPrint();

    int nv;

private:
    CSRMatrix adj;
    aligned_vector<int> degrees;
};
