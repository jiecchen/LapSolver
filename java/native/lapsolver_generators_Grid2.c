#include "lapsolver_generators_Grid2.h"
#include "common.h"

#define getIdx(i,j) (width*(i) + (j))
void doPopulate(int * restrict src, int * restrict dst, double * restrict weight,
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

JNIEXPORT void JNICALL Java_lapsolver_generators_Grid2_populateC
  (JNIEnv *env, jobject cls,
   jintArray jSrc, jintArray jDst, jdoubleArray jWeight,
   jint jHeight, jint jWidth, jint jVerticalWeight)
{
    // convert to c-types
    jint *src = (*env)->GetPrimitiveArrayCritical(env, jSrc, NULL);
    jint *dst = (*env)->GetPrimitiveArrayCritical(env, jDst, NULL);
    jdouble *weight = (*env)->GetPrimitiveArrayCritical(env, jWeight, NULL);

    const int height = (int) jHeight;
    const int width  = (int) jWidth;
    const double verticalWeight = (double) jVerticalWeight;

    doPopulate(src, dst, weight, height, width, verticalWeight);

    (*env)->ReleasePrimitiveArrayCritical(env, jSrc, src, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jDst, dst, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jWeight, weight, 0);
}
