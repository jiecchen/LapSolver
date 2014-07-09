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

    //fills in most of the edges
    #pragma omp parallel for 
    for (int i = 0; i < shortz; i++) {
        for (int j = 0; j < shorty; j++) {
            int e = (i * shorty + j) * (3 * x - 1);
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
        }
    }

    int offset = shortz * shorty * (3 * x - 1);
    //bottom grid (minus the bottom-back edges)
    #pragma omp parallel for
    for (int i = 0; i < shortz; i++) {
        int e = offset + i * (2 * x - 1);
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
    }

    offset += shortz * (2 * x - 1);
    //back grid (minus the bottom-back edges)
    #pragma omp parallel for
    for (int j = 0; j < shorty; j++) {
        int e = offset + j * (2 * x - 1);
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
    }

    int e = offset + shorty * (2 * x - 1);
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

    // the majority of the hypercube
    #pragma omp parallel for 
    for (int i = 0; i < shortz; i++) {
        for (int j = 0; j < shorty; j++) {
            for (int k = 0; k < shortx; k++) {
                int e = (i * shorty * shortx + j * shortx + k) * (4 * w - 1);
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
            }
        }
    }

    int offset = shortz * shorty * shortx * (4 * w - 1);
    // most of i = shortz cube
    #pragma omp parallel for
    for (int j = 0; j < shorty; j++) {
        for (int k = 0; k < shortx; k++) {
            int e = offset + (j * shortx + k) * (3 * w - 1);
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
        }
    }

    offset += shorty * shortx * (3 * w - 1);
    // i = shortz, j = shorty portion of cube
    #pragma omp parallel for
    for (int k = 0; k < shortx; k++) {
        int e = offset + k * (2 * w - 1);
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
    }

    offset += shortx * (2 * w - 1);
    // i = shortz, k = shortx portion of cube
    #pragma omp parallel for
    for (int j = 0; j < shorty; j++) {
        int e = offset + j * (2 * w - 1);
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
    }

    offset += shorty * (2 * w - 1);
    // remainder of i = shortz cube
    #pragma omp parallel for
    for (int l = 0; l < shortw; l++) {
        //w
        src[offset+l] = get4d(shortz, shorty, shortx, l);
        dst[offset+l] = get4d(shortz, shorty, shortx, l + 1);
        weight[offset+l] = 1.0;
    }

    offset += shortw;
    // some of the j = shorty region
    #pragma omp parallel for
    for (int i = 0; i < shortz; i++) {
        for (int k = 0; k < shortx; k++) {
            int e = offset + (i * shortx + k) * (3 * w - 1);
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
        }
    }

    offset += shortz * shortx * (3 * w - 1);
    // remaining j = shorty, k = shortx region
    #pragma omp parallel for 
    for (int i = 0; i < shortz; i++) {
        int e = offset + i * (2 * w - 1);
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
    }

    offset += shortz * (2 * w - 1);
    #pragma omp parallel for 
    for (int i = 0; i < shortz; i++) {
        for (int j = 0; j < shorty; j++) {
            int e = offset + (i * shorty + j) * (3 * w - 1);
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
        }
    }
}

void cbt(int * restrict src, int * restrict dst, double * restrict weight, const int h) {

    int n = (1 << h) - 1;
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        int x = i << 1;
        src[x] = i;
        dst[x] = x + 1;
        weight[x] = 1.0;

        src[x + 1] = i;
        dst[x + 1] = x + 2;
        weight[x + 1] = 1.0;
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
