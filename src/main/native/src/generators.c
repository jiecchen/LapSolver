#include "generators.h"

#define getIdx(i,j) (width*(i) + (j))
void grid2(int * restrict src, int * restrict dst, double * restrict weight,
           const int height, const int width, const double verticalWeight)
{
    const int N = width*height;
    const int shortHeight = height-1;
    const int shortWidth = width-1;

    // populate the majority of edges
    #pragma omp parallel for
    for(int i = 0; i < shortHeight; i++) {
        int e = i * (2 * width - 1);
        for(int j = 0; j < width; j++) {
            src[e+j] = getIdx(i, j);
            dst[e+j] = getIdx(i+1, j);
            weight[e+j] = verticalWeight;
        }
        e += width;
        for(int j = 0; j < shortWidth; j++) {
            src[e+j] = getIdx(i, j);
            dst[e+j] = getIdx(i, j+1);
            weight[e+j] = 1.0;
        }
    }

    // populate bottom edge
    const int base = shortHeight*(2*width-1);
    src += base; dst += base; weight += base;
    for (int j = 0; j < shortWidth; j++) {
        src[j] = getIdx(shortHeight,j);
        dst[j] = getIdx(shortHeight,j+1);
        weight[j] = 1.0;
    }
}

void triangleGrid2(int * restrict src, int * restrict dst, double * restrict weight,
                   const int height, const int width)
{
    grid2(src, dst, weight, height, width, 1.0);

    const int shortHeight = height - 1;
    const int shortWidth = width - 1;
    const int totalEdges = (2 * width * height) - width - height;

    #pragma omp parallel for
    for(int i = 0; i < shortHeight; i++) { // for each row
        int e = totalEdges + i * (width - 1);
        for(int j = 0; j < shortWidth; j++) {
            src[e+j] = getIdx(i,j);
            dst[e+j] = getIdx(i+1,j+1);
            weight[e+j] = 1.0;
        }
    }
}
