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

// private:
    CSRMatrix adj;
};

static void printGraph(const Graph &g) {
    printf("Graph has dimension %d and %d entries\n", g.adj.dim, g.adj.nnz);
    printf("data: ");
    for (double data : g.adj.data)
        printf("%.0f ", data);
    
    printf("\ncols: ");
    for (int cols : g.adj.cols)
        printf("%d ", cols);

    printf("\nrows: ");
    for (int rows : g.adj.rows)
        printf("%d ", rows);

    printf("\n");
}
