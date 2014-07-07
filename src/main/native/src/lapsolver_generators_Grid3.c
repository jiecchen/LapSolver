#include "lapsolver_generators_Grid3.h"
#include "generators.h"

JNIEXPORT void JNICALL Java_lapsolver_generators_Grid3_populateg3
	(JNIEnv *env, jobject cls,
	 jintArray jSrc, jintArray jDst, jdoubleArray jWeight,
	 jint jx, jint jy, jint jz, jint jyweight, jint jzweight)
{
	// convert to c-types
	jint *src = (*env)->GetPrimitiveArrayCritical(env, jSrc, NULL);
	jint *dst = (*env)->GetPrimitiveArrayCritical(env, jDst, NULL);
    jdouble *weight = (*env)->GetPrimitiveArrayCritical(env, jWeight, NULL);

    const int x = (int) jx;
    const int y = (int) jy;
    const int z = (int) jz;
    const double yweight = (double) jyweight;
    const double zweight = (double) jzweight;

    grid3(src, dst, weight, x, y, z, yweight, zweight);

    (*env)->ReleasePrimitiveArrayCritical(env, jSrc, src, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jDst, dst, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jWeight, weight, 0);
}