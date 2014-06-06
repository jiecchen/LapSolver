/**
 * @file KelnerSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An implementation of the primal-dual algorithm of Kelner et al.
 */

package lapsolver.solvers.kelner;

import lapsolver.solvers.Solver;
import lapsolver.lsst.SimulPathTree;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.EdgeList;

public class KelnerSolver implements Solver {
    private Tree spanningTree;
    private PathUpdateTree pathTree;

    private EdgeList offEdges;
    private int[] offLca;
    private double[] offStretch;

    // initialize solver on a particular graph, and perform preprocessing
    @Override
    public void init(Graph graph) {
        // compute LSST
        spanningTree = (new SimulPathTree(graph)).getTree();
    }

    // solve for x in Lx = b
    @Override
    public double[] solve(double[] b) {
        return null;
    }


}
