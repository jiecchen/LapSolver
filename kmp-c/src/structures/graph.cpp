#include <mkl.h>
#include "graph.h"

Graph::Graph(EdgeList &edges)
{
    adj.nnz = edges.ne * 2;
    adj.dim = edges.nv;

    // Temp storage. Will be freed when the constructor exits.
    aligned_vector<double> triData(edges.ne);
    aligned_vector<int> triCols(edges.ne);
    aligned_vector<int> triRows(adj.dim);

    // COO to upper-triangular CSR
    int nnz = edges.ne;
    int info = -1;
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
                &info);

    // Allocate final memory
    aligned_vector<double> data(adj.nnz);
    aligned_vector<int> cols(adj.nnz);
    aligned_vector<int> rows(adj.dim);

    // B = A + A'
    int request = 0, sort = 0;
    double beta = 1.0;
    int WTF_NNZ = adj.nnz + 1; // TODO: wat iz dis, mkl?
    mkl_dcsradd("t", &request, &sort, &adj.dim, &adj.dim,
                triData.data(), triCols.data(), triRows.data(), // A
                &beta, triData.data(), triCols.data(), triRows.data(), // bA'
                data.data(), cols.data(), rows.data(),
                &WTF_NNZ, &info);

    // Convert back to 0-based indexing
    cols.data()[0:adj.nnz]--;
    rows.data()[0:adj.dim]--;

    adj.data = data;
    adj.cols = cols;
    adj.rows = rows;
}
