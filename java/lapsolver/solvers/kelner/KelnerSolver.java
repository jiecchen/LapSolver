/**
 * @file KelnerSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An implementation of the primal-dual algorithm of Kelner et al.
 */

package lapsolver.solvers.kelner;

import lapsolver.algorithms.DiscreteSampler;
import lapsolver.algorithms.Stretch;
import lapsolver.solvers.Solver;
import lapsolver.lsst.SimulPathTree;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.EdgeList;
import lapsolver.lsst.SpanningTreeStrategy;
import lapsolver.util.TreeUtils;

public class KelnerSolver implements Solver {
    // graph-theoretic data structures
    private Graph graph;
    private Tree spanningTree;
    private SpanningTreeStrategy treeStrategy;
    private PathUpdateTree pathTree;

    // edge data to be preprocessed
    private EdgeList offEdges;
    private int[] offLca;
    private double[] offStretch;
    private DiscreteSampler edgeSampler;

    // initialize solver with a spanning tree strategy
    public KelnerSolver(SpanningTreeStrategy treeStrategy) {
        this.treeStrategy = treeStrategy;
    }

    // initialize solver on a particular graph, and perform preprocessing
    @Override
    public void init(Graph graph) {
        // compute LSST
        this.graph = graph;
        spanningTree = treeStrategy.getTree(graph);

        // get off-tree edges; compute their stretches
        offEdges = TreeUtils.getOffTreeEdges(graph, spanningTree);
        offStretch = Stretch.compute(graph, spanningTree, offEdges).allStretches;

        


    }

    // solve for x in Lx = b
    @Override
    public double[] solve(double[] b) {
        return null;
    }


}
