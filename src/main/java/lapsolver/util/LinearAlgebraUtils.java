/**
 * @file LinearAlgebraUtils.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Mon Jul 21 2014
 *
 * Static methods for common linear algebra operations.
 */

package lapsolver.util;

public class LinearAlgebraUtils {
    /**
     * Computes the dot product between two vectors.
     */
    public static double dot (double[] u, double[] v) {
        double result = 0;
        for (int i = 0; i < u.length; i++) {
            result += u[i] * v[i];
        }
        return result;
    }

    /**
     * Computes the norm of a vector.
     */
    public static double norm (double[] u) {
        double result = 0;
        for (int i = 0; i < u.length; i++) {
            result += u[i] * u[i];
        }
        return Math.sqrt(result);
    }

    /**
     * Normalize a vector in place.
     */
    public static void normalize (double[] v) {
        double norm = norm(v);
        for (int i = 0; i < v.length; i++) {
            v[i] /= norm;
        }
    }

    /**
     * Computes the sum of two arrays, element wise
     */
     public static double[] add (double[] u, double[] v) {
         double[] res = new double[u.length];
         for (int i = 0; i < res.length; i++)
             res[i] = u[i] + v[i];
         return res;
     }

    /**
     * Computes the sum of two arrays, element wise
     */
     public static double[] subtract (double[] u, double[] v) {
         double[] res = new double[u.length];
         for (int i = 0; i < res.length; i++)
             res[i] = u[i] - v[i];
         return res;
     }

    /**
     * Computes the sum of two arrays, element wise
     */
     public static double[] scale (double[] u, double x) {
         double[] res = new double[u.length];
         for (int i = 0; i < res.length; i++)
             res[i] = u[i] * x;
         return res;
     }
}
