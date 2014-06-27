/**
 * @file Solver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An interface for Laplacian solvers.
 */

package lapsolver.solvers;

import lapsolver.Graph;

public interface Solver {
    /**
     * Initialize solver on a particular graph, and perform preprocessing.
     * @param graph The input graph (weights are conductances).
     */
    public void init(Graph graph);

    // solve for x in Lx = b
    public double[] solve(double[] b);
}
