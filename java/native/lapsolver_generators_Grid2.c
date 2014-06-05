#include "lapsolver_generators_Grid2.h"
#include <string.h>
#define getIdx(i,j) (width*(i) + (j))

JNIEXPORT void JNICALL Java_lapsolver_generators_Grid2_populateC(
   JNIEnv *env, jobject cls,
   jintArray jSrc, jintArray jDst, jdoubleArray jWeight, 
   jint jHeight, jint jWidth, jint jVerticalWeight)
{
    // convert to c-types
    const int height = (int) jHeight;
    const int width  = (int) jWidth;
    const int verticalWeight = (int) jVerticalWeight;

    int *src = (*env)->GetIntArrayElements(env, jSrc, 0);
    int *dst = (*env)->GetIntArrayElements(env, jDst, 0);
    double *weight = (*env)->GetDoubleArrayElements(env, jWeight, 0);

    const int N = width*height;

    int e = 0;
    const int shortHeight = height-1;
    const int shortWidth = width-1;

    // populate the majority of edges
    for(int i = 0; i < shortHeight; i++) {
        for(int j = 0; j < shortWidth; j++, e++) { // vertical - vectorized
            src[e] = getIdx(i, j);
            dst[e] = getIdx(i+1, j);
            weight[e] = verticalWeight;
        }
        for(int j = 0; j < shortWidth; j++, e++) { // horizontal - vectorized
            src[e] = getIdx(i, j);
            dst[e] = getIdx(i, j+1);
            weight[e] = 1.0;
        }
    }
    // populate right edge - not vectorized
    const int edgeBound = N-1;
    for (int i = width-1; i < edgeBound; i+=width, e++) {
        src[e] = i;
        dst[e] = i+width;
        weight[e] = verticalWeight;
    }

    // populate bottom edge - vectorized
    for (int j = width*(height-1); j < edgeBound; j++, e++) {
        src[e] = j;
        dst[e] = j+1;
        weight[e] = 1.0;
    }

    (*env)->ReleaseIntArrayElements(env, jSrc, src, 0);
    (*env)->ReleaseIntArrayElements(env, jDst, dst, 0);
    (*env)->ReleaseDoubleArrayElements(env, jWeight, weight, 0);
}
