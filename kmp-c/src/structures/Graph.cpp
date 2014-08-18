#include <utility>
#include <mkl.h>
#include <cilk/cilk.h>
#include <math.h>
#include "Graph.h"

Graph::Graph() { }

Graph::Graph(const Graph &rval)
    : adj(rval.adj),
      degrees(rval.degrees)
{

}

Graph::Graph(Graph &&rval)
    : adj(std::move(rval.adj)),
      degrees(std::move(rval.degrees))
{

}

Graph &Graph::operator=(const Graph &g)
{
    if (this != &g)
    {
        adj = g.adj;
        degrees = g.degrees;
    }
    return *this;
}

Graph &Graph::operator=(Graph && g)
{
    adj = std::move(g.adj);
    degrees = std::move(g.degrees);
    return *this;
}

#define mkl_throw(f) throw fprintf(stderr, "Error! " f " failed with code %d\n", error), error;

Graph::Graph(EdgeList &&edges)
{
    adj.nnz = edges.ne * 2;
    adj.dim = edges.nv;

    // Temp storage. Will be freed when the constructor exits.
    aligned_vector<double> triData(edges.ne);
    aligned_vector<int> triCols(edges.ne);
    aligned_vector<int> triRows(adj.dim + 1);

    // COO to upper-triangular CSR
    int nnz = edges.ne;
    int error = 0;
    int job [] =
    {
        2,   // COO is converted to the CSR format; the column indices are increasing within each row
        1,   // one-based indexing for the matrix in CSR format is used (required by DSCRADD)
        0,   // zero-based indexing for the matrix in coordinate format is used
        0,   // ???
        nnz, // number of the non-zero elements of the matrix A
        0    // all arrays acsr, ja, ia are filled in for the output storage
    };
    mkl_dcsrcoo(job, &adj.dim,
                triData.data(),
                triCols.data(),
                triRows.data(),
                &nnz,
                edges.w.data(),
                edges.u.data(),
                edges.v.data(),
                &error);

    if (error)
        mkl_throw("mkl_dcsrcoo");

    // Allocate final memory
    aligned_vector<int> rows(adj.dim + 1);
    aligned_vector<double> data(adj.nnz);
    aligned_vector<int> cols(adj.nnz);

    // B = A + A'
    int request = 1, sort = 0;
    double beta = 1.0;
    mkl_dcsradd("t", &request, &sort, &adj.dim, &adj.dim,
                triData.data(), triCols.data(), triRows.data(), // A
                &beta, triData.data(), triCols.data(), triRows.data(), // bA'
                NULL, NULL, rows.data(),
                NULL, &error);

    request = 2;
    mkl_dcsradd("t", &request, &sort, &adj.dim, &adj.dim,
                triData.data(), triCols.data(), triRows.data(), // A
                &beta, triData.data(), triCols.data(), triRows.data(), // bA'
                data.data(), cols.data(), rows.data(),
                NULL, &error);

    if (error)
        mkl_throw("mkl_dcsradd");

    // Convert back to 0-based indexing
    cols.data()[0:adj.nnz]--;
    rows.data()[0:adj.dim + 1]--;

    adj.data = data;
    adj.cols = cols;
    adj.rows = rows;

    degrees = aligned_vector<int>(adj.dim);
    degrees.data()[0:adj.dim] = rows.data()[1:adj.dim] - rows.data()[0:adj.dim];
}

void Graph::debugPrint()
{
    printf("Graph has dimension %d and %d entries\n", adj.dim, adj.nnz);
    printf("data: ");
    for (double data : adj.data)
        printf("%.0f ", data);

    printf("\ncols: ");
    for (int cols : adj.cols)
        printf("%d ", cols);

    printf("\nrows: ");
    for (int rows : adj.rows)
        printf("%d ", rows);

    printf("\ndegs: ");
    for (int deg : degrees)
        printf("%d ", deg);

    printf("\n");
}
