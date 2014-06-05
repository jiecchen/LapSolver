/**
 * @file KelnerSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An implementation of the primal-dual algorithm of Kelner et al.
 */

package lapsolver.solvers.kelner;

import lapsolver.solvers.Solver;
import lapsolver.lsst.SpanningTreeStrategy;
import lapsolver.Graph;
import lapsolver.Tree;

public class KelnerSolver implements Solver {
    private SpanningTreeStrategy spanningTreeStrategy;
    private Tree spanningTree;
    private PathUpdateTree pathTree;

    // initialize solver on a particular graph, and perform preprocessing
    @Override
    public void init(Graph g) {

    }

    // solve for x in Lx = b
    @Override
    public double[] solve(double[] b) {
        return null;
    }
}
