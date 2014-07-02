/**
 * @file Solver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An interface for Laplacian solvers.
 */

package lapsolver.solvers;

import lapsolver.Graph;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.apache.commons.math3.linear.RealVector;

public abstract class Solver extends RealLinearOperator {
    public Graph graph;

    /**
     * Initialize solver on a particular graph, and perform preprocessing.
     * @param graph The input graph (weights are conductances).
     * @param d Excess diagonal entries (null if we want to solve a Laplacian system).
     */
    public abstract void init (Graph graph, double[] d);

    /**
     * Initialize solver on a particular graph, and perform preprocessing. Pass in null for excess diagonal entries.
     * @param graph The input graph (weights are conductances).
     */
    public void init (Graph graph) {
        init (graph, null);
    }

    /**
     * Solve the system (L+D)x = b.
     * @param b The boundary condition b.
     * @return The solution x.
     */
    public abstract double[] solve (double[] b);

    @Override
    public int getRowDimension() {
        return graph.nv;
    }

    @Override
    public int getColumnDimension() {
        return graph.nv;
    }

    @Override
    public RealVector operate(RealVector x) throws DimensionMismatchException {
        return new ArrayRealVector( solve(x.toArray()) );
    }
}
