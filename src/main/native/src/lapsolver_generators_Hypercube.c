#include "lapsolver_generators_Hypercube.h"
#include "generators.h"

JNIEXPORT void JNICALL Java_lapsolver_generators_Hypercube_populateg4
	(JNIEnv *env, jobject cls,
	 jintArray jSrc, jintArray jDst, jdoubleArray jWeight, jint jw,
	 jint jx, jint jy, jint jz, jdouble jxweight, jdouble jyweight, jdouble jzweight)
{
	// convert to c-types
	jint *src = (*env)->GetPrimitiveArrayCritical(env, jSrc, NULL);
	jint *dst = (*env)->GetPrimitiveArrayCritical(env, jDst, NULL);
    jdouble *weight = (*env)->GetPrimitiveArrayCritical(env, jWeight, NULL);

    const int w = (int) jw;
    const int x = (int) jx;
    const int y = (int) jy;
    const int z = (int) jz;
    const double xweight = (double) jxweight;
    const double yweight = (double) jyweight;
    const double zweight = (double) jzweight;

    hypercube(src, dst, weight, w, x, y, z, xweight, yweight, zweight);

    (*env)->ReleasePrimitiveArrayCritical(env, jSrc, src, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jDst, dst, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, jWeight, weight, 0);
}




