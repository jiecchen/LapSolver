#include "lapsolver_generators_cbt.h"
#include "generators.h"

JNIEXPORT void JNICALL Java_lapsolver_generators_CBT_populateCBT
  (JNIEnv *env, jobject cls,
   jintArray jSrc, jintArray jDst, jdoubleArray jWeight, jint jh)
{
    // convert to c-types
    jint *src = (*env)->GetPrimitiveArrayCritical(env, jSrc, NULL);
    jint *dst = (*env)->GetPrimitiveArrayCritical(env, jDst, NULL);
    jdouble *weight = (*env)->GetPrimitiveArrayCritical(env, jWeight, NULL);

    const int h = (int) jh;

    cbt(src, dst, weight, h);

    (*env)->ReleasePrimitiveArrayCritical(env, jSrc, src, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jDst, dst, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jWeight, weight, 0);
}