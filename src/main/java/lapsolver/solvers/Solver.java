/**
 * @file Solver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An interface for Laplacian solvers.
 */

package lapsolver.solvers;

import lapsolver.Graph;
import lapsolver.LinearOperator;

public abstract class Solver implements LinearOperator {
    public Graph graph;
    public double[] d;

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

    public Solver initialize(Graph graph, double[] d) {
        init(graph, d);
        return this;
    }

    public Solver initialize(Graph graph) {
        return initialize(graph, null);
    }

    /**
     * Solve the system (L+D)x = b.
     * @param b The boundary condition b.
     * @return The solution x.
     */
    public abstract double[] solve (double[] b);

    @Override
    public double[] apply(double[] b) {
        return solve(b);
    }
}
