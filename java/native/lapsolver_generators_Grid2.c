#include "lapsolver_generators_Grid2.h"
#include <string.h>
#define getIdx(i,j) (width*(i) + (j))

JNIEXPORT void JNICALL Java_lapsolver_generators_Grid2_populateC(
   JNIEnv *env, jobject cls,
   jobject jSrc, jobject jDst, jobject jWeight,
   jint jHeight, jint jWidth, jint jVerticalWeight)
{
    // convert to c-types
    jint *src = (*env)->GetDirectBufferAddress(env, jSrc);
    jint *dst = (*env)->GetDirectBufferAddress(env, jDst);
    jdouble *weight = (*env)->GetDirectBufferAddress(env, jWeight);

    if(!src || !dst || !weight)
        return;

    const int height = (int) jHeight;
    const int width  = (int) jWidth;
    const double verticalWeight = (double) jVerticalWeight;

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
}

