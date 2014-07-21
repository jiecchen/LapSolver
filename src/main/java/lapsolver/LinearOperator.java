/**
 * @file LinearOperator.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Mon Jul 21 2014
 *
 * An interface for linear operators.
 */

package lapsolver;

public interface LinearOperator {
    public double[] apply(double[] x);
}
