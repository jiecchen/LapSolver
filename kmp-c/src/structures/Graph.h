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
    Graph(const Graph &rval);
    Graph(Graph &&rval);

    Graph &operator=(const Graph &g);
    Graph &operator=(Graph && g);

    explicit Graph(EdgeList &&edges);

    inline double& weight(int v, int i)
    {
        return adj.data[adj.rows[v] + i];
    }

    inline int& neighbor(int v, int i)
    {
        return adj.cols[adj.rows[v] + i];
    }

    inline const double& weight(int v, int i) const
    {
        return adj.data[adj.rows[v] + i];
    }

    inline const int& neighbor(int v, int i) const
    {
        return adj.cols[adj.rows[v] + i];
    }

    inline int degree(int v) const
    {
        return degrees[v];
    }

    inline int dataIndex(int v, int i) const
    {
        return adj.rows[v] + i;
    }

    void debugPrint();

    inline int nv() const
    {
        return adj.dim;
    }

    inline int ne() const
    {
        return adj.nnz >> 1;
    }

private:
    CSRMatrix<RealGeneral> adj;
    aligned_vector<int> degrees;
};
