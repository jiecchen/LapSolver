#include <iostream>
#include <vector>
#include <cstdlib>
#include <sys/time.h>
#include <ctime>
#include <numeric>
#include <mkl.h>
#include <util/aligned.h>
#include <cmath>
#include <util/Benchmark.h>
using namespace std;

#define RANDOM (((double) rand()) / RAND_MAX)

int main(int argc, char *argv[])
{
    vmlSetMode(VML_EP);
    atexit(MKL_Free_Buffers);
    srand(time(NULL));

    const int SIZE = 1000;
    aligned_vector<double> x(SIZE);
    aligned_vector<double> y(SIZE);
    for (int i = 0; i < SIZE; ++i)
    {
        x[i] = RANDOM;
        y[i] = RANDOM;
    }

    aligned_vector<double> z;
    auto bench = make_benchmark(argc, argv, [&] () {
        z = x + y;
    });
    z.data()[0:SIZE] -= x.data()[0:SIZE];
    z.data()[0:SIZE] -= y.data()[0:SIZE];
    for (auto &pt : z)
        printf("%.2f\n", fabs(pt));

    return 0;
}
