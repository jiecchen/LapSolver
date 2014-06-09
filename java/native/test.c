#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

int main()
{
    const int height = 1000;
    const int width  = 1000;
    const int ne = (2 * width * height) - width - height;
    int *src = aligned_alloc(32, sizeof(int)*ne);
    int *dst = aligned_alloc(32, sizeof(int)*ne);
    double *weights = aligned_alloc(32, sizeof(double)*ne);

    doPopulate(src, dst, weights, height, width, 1.0);

    int ret = src[ne-1] + dst[ne-1];
    free(src);
    free(dst);
    free(weights);
    return ret;
}
