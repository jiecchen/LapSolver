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

void grid3(int * restrict src, int * restrict dst, double * restrict weight,
           const int x, const int y, const int z,
           const double yweight, const double zweight)
{
    const int shortx = x - 1;
    const int shorty = y - 1;
    const int shortz = z - 1;
    const int slice = x * shorty + shortx * shorty;

    int e = 0; //edge number counter

    //fills in most of the edges
    #pragma omp parallel for 
    for (int i = 0; i < shortz; i++) {
        for (int j = 0; j < shorty; j++) {
            //edges in y and z directions
            for (int k = 0; k < x; k++) {
                //y
                src[e+k] = i * x * y + j * x + k;
                dst[e+k] = i * x * y + (j + 1) * y + k;
                weight[e+k] = yweight;
                //z
                src[e+k+1] = i * x * y + j * x + k;
                dst[e+k+1] = (i + 1) * x * y + j * x + k;
                weight[e+k+1] = zweight;
                e++;
            }
            e += x;
            //edges in x direction
            for (int k = 0; k < shortx; k++) {
                src[e+k] = i * x * y + j * x + k;
                dst[e+k] = i * x * y + j * x + k + 1;
                weight[e+k] = 1.0;
            }
            e += shortx;
        }
    }
    //bottom grid (minus the bottom-back edges)
    #pragma omp parallel for
    for (int i = 0; i < shortz; i++) {
        //edges in z direction
        for (int k = 0; k < x; k++) {
            src[e+k] = i * x * y + shorty * x + k;
            dst[e+k] = (i + 1) * x * y + shorty * x + k;
            weight[e+k] = zweight;
        }
        e += x;
        //edges in x direction
        for (int k = 0; k < shortx; k++) {
            src[e+k] = i * x * y + shorty * x + k;
            dst[e+k] = i * x * y + shorty * x + k + 1;
            weight[e+k] = 1.0;
        }
        e += shortx;
    }
    //back grid (minus the bottom-back edges)
    #pragma omp parallel for
    for (int j = 0; j < shorty; j++) {
        //y direction
        for (int k = 0; k < x; k++) {
            src[e+k] = shortz * x * y + j * x + k;
            dst[e+k] = shortz * x * y + (j + 1) * y + k;
            weight[e+k] = yweight;
        }
        e += x;
        //x direction
        for (int k = 0; k < shortx; k++) {
            src[e+k] = shortz * x * y + j * x + k;
            dst[e+k] = shortz * x * y + j * x + k + 1;
            weight[e+k] = 1.0;
        }
        e += shortx;
    }
    //bottom-back edges
    #pragma omp parallel for
    for (int k = 0; k < shortx; k++) {
        src[e+k] = shortz * x * y + shorty * x + k;
        dst[e+k] = shortz * x * y + shorty * x + k + 1;
        weight[e+k] = 1.0;
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
