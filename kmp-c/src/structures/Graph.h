/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Graph class for easy use with BLAS
 */
#pragma once
#include <utility>
#include "util/aligned.h"
#include "EdgeList.h"
#include "CSRMatrix.h"

class Graph
{
public:
    Graph();
    Graph(const Graph &rval);
    Graph(Graph &&rval);

    Graph &operator=(const Graph &g);
    Graph &operator=(Graph && g);

    explicit Graph(EdgeList &&edges);

    inline const double *getWeights(int v) const
    {
        return adj.data.data() + adj.rows[v];
    }

    inline const int *getNeighbors(int v) const
    {
        return adj.cols.data() + adj.rows[v];
    }

    inline int getDegree(int v) const
    {
        return degrees[v];
    }

    inline int getDataIndex(int v, int i) const
    {
        return adj.rows[v] + i;
    }

    void debugPrint();

    int nv;

private:
    CSRMatrix adj;
    aligned_vector<int> degrees;
};
