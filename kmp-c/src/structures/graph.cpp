#include <mkl.h>
#include "graph.h"

Graph::Graph(EdgeList &edges)
{
    adj.nnz = edges.ne * 2;
    adj.dim = edges.nv;
    aligned_vector<double> triData(edges.ne);
    aligned_vector<int> triCols(edges.ne);
    aligned_vector<int> triRows(adj.dim);

    int dim = edges.nv;
    int nnz = edges.ne;
    int info = -1;
    int job [] =
    {
        2,   // COO is converted to the CSR format; the column indices are increasing within each row
        1,   // zero-based indexing for the matrix in CSR format is used
        0,   // zero-based indexing for the matrix in coordinate format is used
        0,   // ???
        nnz, // number of the non-zero elements of the matrix A
        0    // all arrays acsr, ja, ia are filled in for the output storage
    };
    mkl_dcsrcoo(job,
                &dim,
                triData.data(),
                triCols.data(),
                triRows.data(),
                &nnz,
                edges.w.data(),
                edges.u.data(),
                edges.v.data(),
                &info);

    aligned_vector<double> data(adj.nnz);
    aligned_vector<int> cols(adj.nnz);
    aligned_vector<int> rows(adj.dim);

    int WTF_NNZ = adj.nnz + 1; // TODO: wat iz dis, mkl?
    mkl_dcsradd("t", &(int){0}, &(int){0}, &dim, &dim,
                triData.data(), triCols.data(), triRows.data(), // A
                &(double){1.0},
                triData.data(), triCols.data(), triRows.data(), // bA'
                data.data(), cols.data(), rows.data(),
                &WTF_NNZ, &info);

    adj.data = data;
    adj.cols = cols;
    adj.rows = rows;
}
