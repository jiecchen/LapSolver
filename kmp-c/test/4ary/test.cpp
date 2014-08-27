#include <vector>
#include <cstdio>
#include "structures/heap/Aligned4aryHeap.h"
#include "util/Benchmark.h"
#include "mkl.h"
using namespace std;

int main(int argc, char **argv)
{
    int nElements, v; double k;
    scanf("%d\n", &nElements);

    vector<double> loadedKeys(nElements);
    vector<int> loadedValues(nElements);
    vector<int> sortedValues(nElements);

    for (int i = 0; i < nElements; ++i)
    {
        scanf("%lf %d\n", &k, &v);
        loadedKeys[i] = k;
        loadedValues[i] = v;
    }

    auto bench = make_benchmark(argc, argv, [&] () {
        Aligned4aryHeap<double, int> heap(nElements);
        for(unsigned i = 0; i < nElements; ++i)
            heap.push(loadedKeys[i], loadedValues[i]);
        for(unsigned i = 0; i < nElements; ++i)
            heap.pop(&k, &sortedValues[i]);
    });

    for(unsigned i = 0; i < nElements; ++i)
        printf("%d\n", sortedValues[i]);

    mkl_free_buffers();
    return 0;
}
