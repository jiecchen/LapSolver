#include "lapsolver_generators_Grid2.h"

#define getIdx(i,j) (width*(i) + (j))
void doPopulate(int * restrict src, int * restrict dst, double * restrict weight,
                const int height, const int width, const double verticalWeight)
{
    const int N = width*height;
    const int shortHeight = height-1;
    const int shortWidth = width-1;

    int e = 0;
    // populate the majority of edges
    for(int i = 0; i < shortHeight; i++) {
        for(int j = 0; j < width; j++, e++) {
            src[e] = getIdx(i, j);
            dst[e] = getIdx(i+1, j);
            weight[e] = verticalWeight;
        }
        for(int j = 0; j < shortWidth; j++, e++) {
            src[e] = getIdx(i, j);
            dst[e] = getIdx(i, j+1);
            weight[e] = 1.0;
        }
    }

    // populate bottom edge
    for (int j = width*shortHeight+1; j < N; j++, e++) {
        src[e] = j-1;
        dst[e] = j;
        weight[e] = 1.0;
    }
}

//JNIEXPORT void JNICALL Java_lapsolver_generators_Grid2_populateC(
//   JNIEnv *env, jobject cls,
//   jobject jSrc, jobject jDst, jobject jWeight,
//   jint jHeight, jint jWidth, jint jVerticalWeight)
//{
//    // convert to c-types
//    jint *src = (*env)->GetDirectBufferAddress(env, jSrc);
//    jint *dst = (*env)->GetDirectBufferAddress(env, jDst);
//    jdouble *weight = (*env)->GetDirectBufferAddress(env, jWeight);
//
//    if(!src || !dst || !weight)
//        return;
//
//    const int height = (int) jHeight;
//    const int width  = (int) jWidth;
//    const double verticalWeight = (double) jVerticalWeight;
//
//    doPopulate(src, dst, weight, height, width, verticalWeight);
//}

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
