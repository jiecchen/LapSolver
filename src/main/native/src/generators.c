#include "generators.h"

#define getIdx(i,j) (width*(i) + (j))
#define get4d(i, j, k, l) (i * w * x * y + j * w * x + k * w + l) //from 4d coordinates to 1d

void grid2(int * restrict src, int * restrict dst, double * restrict weight,
           const int height, const int width, const double verticalWeight)
{

    const int shortHeight = height-1;
    const int shortWidth = width-1;

    // populate the majority of edges
    #pragma omp parallel for
    for(int i = 0; i < shortHeight; i++) {
        int e = i * (2 * width - 1);
        //y
        for(int j = 0; j < width; j++) {
            src[e+j] = getIdx(i, j);
            dst[e+j] = getIdx(i+1, j);
            weight[e+j] = verticalWeight;
        }
        e += width;
        //x
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
                e++;
                //z
                src[e+k] = i * x * y + j * x + k;
                dst[e+k] = (i + 1) * x * y + j * x + k;
                weight[e+k] = zweight;
                
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

void hypercube(int * restrict src, int * restrict dst, double * restrict weight,
               const int w, const int x, const int y, const int z,
               const double xweight, const double yweight, const double zweight) {

    const int shortw = w - 1;
    const int shortx = x - 1;
    const int shorty = y - 1;
    const int shortz = z - 1;

    int e = 0;

    // the majority of the hypercube
    #pragma omp parallel for 
    for (int i = 0; i < shortz; i++) {
        for (int j = 0; j < shorty; j++) {
            for (int k = 0; k < shortx; k++) {
                for (int l = 0; l < w; l++) {
                    //z
                    src[e+l] = get4d(i, j, k, l);
                    dst[e+l] = get4d(i + 1, j, k, l);
                    weight[e+l] = zweight;
                    e++;
                    //y
                    src[e+l] = get4d(i, j, k, l);
                    dst[e+l] = get4d(i, j + 1, k, l);
                    weight[e+l] = yweight;
                    e++;
                    //x
                    src[e+l] = get4d(i, j, k, l);
                    dst[e+l] = get4d(i, j, k + 1, l);
                    weight[e+l] = xweight;
                }
                e += w;
                for (int l = 0; l < shortw; l++) {
                    //w
                    src[e+l] = get4d(i, j, k, l);
                    dst[e+l] = get4d(i, j, k, l + 1);
                    weight[e+l] = 1.0;
                }
                e += shortw;
            }
        }
    }

    // most of i = shortz cube
    #pragma omp parallel for
    for (int j = 0; j < shorty; j++) {
        for (int k = 0; k < shortx; k++) {
            for (int l = 0; l < w; l++) {
                //y
                src[e+l] = get4d(shortz, j, k, l);
                dst[e+l] = get4d(shortz, j + 1, k, l);
                weight[e+l] = yweight;
                e++;
                //x
                src[e+l] = get4d(shortz, j, k, l);
                dst[e+l] = get4d(shortz, j, k + 1, l);
                weight[e+l] = xweight;
            }
            e += w;
            for (int l = 0; l < shortw; l++) {
                //w
                src[e+l] = get4d(shortz, j, k, l);
                dst[e+l] = get4d(shortz, j, k, l + 1);
                weight[e+l] = 1.0;
            }
            e += shortw;
        }
    }

    // i = shortz, j = shorty portion of cube
    #pragma omp parallel for
    for (int k = 0; k < shortx; k++) {
        for (int l = 0; l < w; l++) {
            //x
            src[e+l] = get4d(shortz, shorty, k, l);
            dst[e+l] = get4d(shortz, shorty, k + 1, l);
            weight[e+l] = xweight;
        }
        e += w;
        for (int l = 0; l < shortw; l++) {
            //w
            src[e+l] = get4d(shortz, shorty, k, l);
            dst[e+l] = get4d(shortz, shorty, k, l + 1);
            weight[e+l] = 1.0;
        }
        e += shortw;
    }

    // i = shortz, k = shortx portion of cube
    #pragma omp parallel for
    for (int j = 0; j < shorty; j++) {
        for (int l = 0; l < w; l++) {
            //y
            src[e+l] = get4d(shortz, j, shortx, l);
            dst[e+l] = get4d(shortz, j + 1, shortx, l);
            weight[e+l] = yweight;
        }
        e += w;
        for (int l = 0; l < shortw; l++) {
            //w
            src[e+l] = get4d(shortz, j, shortx, l);
            dst[e+l] = get4d(shortz, j, shortx, l + 1);
            weight[e+l] = 1.0;
        }
        e += shortw;
    }

    // remainder of i = shortz cube
    #pragma omp parallel for
    for (int l = 0; l < shortw; l++) {
        //w
        src[e+l] = get4d(shortz, shorty, shortx, l);
        dst[e+l] = get4d(shortz, shorty, shortx, l + 1);
        weight[e+l] = 1.0;
    }
    e += shortw;

    // some of the j = shorty region
    #pragma omp parallel for
    for (int i = 0; i < shortz; i++) {
        for (int k = 0; k < shortx; k++) {
            for (int l = 0; l < w; l++) {
                //z
                src[e+l] = get4d(i, shorty, k, l);
                dst[e+l] = get4d(i + 1, shorty, k, l);
                weight[e+l] = zweight;
                e++;
                //x
                src[e+l] = get4d(i, shorty, k, l);
                dst[e+l] = get4d(i, shorty, k + 1, l);
                weight[e+l] = xweight;
            }
            e += w;
            for (int l = 0; l < shortw; l++) {
                //w
                src[e+l] = get4d(i, shorty, k, l);
                dst[e+l] = get4d(i, shorty, k, l + 1);
                weight[e+l] = 1.0;
            }
            e += shortw;
        }
    }

    // remaining j = shorty, k = shortx region
    #pragma omp parallel for 
    for (int i = 0; i < shortz; i++) {erm
        for (int l = 0; l < w; l++) {
            //z
            src[e+l] = get4d(i, shorty, shortx, l);
            dst[e+l] = get4d(i + 1, shorty, shortx, l);
            weight[e+l] = 1.0;
        }
        e += w;
        for (int l = 0; l < shortw; l++) {
            //w
            src[e+l] = get4d(i, shorty, shortx, l);
            dst[e+l] = get4d(i, shorty, shortx, l + 1);
            weight[e+l] = 1.0;
        }
        e += shortw;
    }

    #pragma omp parallel for 
    for (int i = 0; i < shortz; i++) {
        for (int j = 0; j < shorty; j++) {
            for (int l = 0; l < w; l++) {
                //z
                src[e+l] = get4d(i, j, shortx, l);
                dst[e+l] = get4d(i + 1, j, shortx, l);
                weight[e+l] = zweight;
                e++;
                //y
                src[e+l] = get4d(i, j, shortx, l);
                dst[e+l] = get4d(i, j + 1, shortx, l);
                weight[e+l] = yweight;
            }
            e += w;
            for (int l = 0; l < shortw; l++) {
                //w
                src[e+l] = get4d(i, j, shortx, l);
                dst[e+l] = get4d(i, j, shortx, l + 1);
                weight[e+l] = 1.0;
            }
            e += shortw;
        }
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
