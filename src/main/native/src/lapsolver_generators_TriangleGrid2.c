#include "lapsolver_generators_Grid2.h"
#include "generators.h"

JNIEXPORT void JNICALL Java_lapsolver_generators_TriangleGrid2_populateC
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

    triangleGrid2(src, dst, weight, height, width);

    (*env)->ReleasePrimitiveArrayCritical(env, jSrc, src, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jDst, dst, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jWeight, weight, 0);
}
