// Test for Cholesky elimination order

#include <cstdlib>
#include <string>
#include <memory>
#include <structures/Graph.h>
#include <structures/GraphLoader.h>
#include <algorithms/PartialCholeskyFactorization.h>
#include <util/Benchmark.h>

int size = 0;

void write(double mat[size][size]) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            printf("%.5lf ", mat[i][j]);
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    atexit(mkl_free_buffers);
    Graph g = GraphLoader::fromStdin();

    aligned_vector <double> diag_values;
    for (int i = 0; i < g.nv(); i++) {
        double x;
        scanf("%lf", &x);
        diag_values.push_back(x);
    }

    int num_steps;
    scanf("%d", &num_steps);

    double L[g.nv()][g.nv()];
    double D[g.nv()][g.nv()];
    memset(L, 0, sizeof(L));
    memset(D, 0, sizeof(D));

    for (int i = 0; i < g.nv(); i++)
        for (int j = 0; j < g.nv(); j++)
            scanf("%lf", &L[i][j]);

    for (int i = 0; i < g.nv(); i++)
        for (int j = 0; j < g.nv(); j++)
            scanf("%lf", &D[i][j]);

    std::shared_ptr<PartialCholeskyFactorization> ldl;

    auto bench = make_benchmark(argc, argv, [&] () {
         ldl.reset(new PartialCholeskyFactorization(g, diag_values, num_steps));
    });

    double myL[g.nv()][g.nv()];
    double myD[g.nv()][g.nv()];
    memset(myL, 0, sizeof(myL));
    memset(myD, 0, sizeof(myD));

    for (int i = 0; i < ldl->L->ne; i++)
        myL[ldl->L->u[i]][ldl->L->v[i]] += ldl->L->w[i];

    for (int i = 0; i < ldl->D->ne; i++)
        myD[ldl->D->u[i]][ldl->D->v[i]] += ldl->D->w[i];

    int eps = 1e-5;
    bool ok = 1;
    for (int i = 0; i < g.nv(); i++)
        for (int j = 0; j < g.nv(); j++) {
            if (abs(L[i][j] - myL[i][j]) > eps || abs(D[i][j] - myD[i][j]) > eps) {
                ok = 0;
                break;
            }
        }

    if (ok)
        printf("OK\n");
    else {
        size = g.nv();
        printf("Matrices differ!\n\n");
        printf("Original L:\n");
        write(L);
        printf("\nOriginal D:\n");
        write(D);
        printf("\nMY L:\n");
        write(myL);
        printf("\nMY D:\n");
        write(myD);
    }

    return 0;
}
