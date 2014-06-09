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
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.EdgeList;
import lapsolver.lsst.SpanningTreeStrategy;
import lapsolver.util.TreeUtils;
import lapsolver.solvers.TreeSolver;

public class KelnerSolver implements Solver {
    // graph data structures
    private Graph graph;
    private Tree spanningTree;
    private SpanningTreeStrategy treeStrategy;
    private PathUpdateTree pathTree;
    private TreeSolver treeSolver;

    // edge data to be preprocessed
    private EdgeList offEdges;
    private int[] offLca;
    private double[] offStretch;
    private DiscreteSampler edgeSampler;

    // algorithm state
    private double[] currentB;
    private double[] currentFlow;
    private int[] order;

    // initialize solver with a spanning tree strategy
    public KelnerSolver(SpanningTreeStrategy treeStrategy) {
        this.treeStrategy = treeStrategy;
    }

    // initialize solver on a particular graph, and perform preprocessing
    @Override
    public void init(Graph graph) {
        // compute LSST, cache BFS order
        this.graph = graph;
        spanningTree = treeStrategy.getTree(graph);
        order = TreeUtils.bfsOrder(spanningTree);

        // get off-tree edges, find stretches, initialize sampler
        offEdges = TreeUtils.getOffTreeEdges(graph, spanningTree);
        offStretch = Stretch.compute(graph, spanningTree, offEdges).allStretches;
        edgeSampler = new DiscreteSampler(offStretch);

        // initialize feasible flow finder for spanning tree
        treeSolver = new TreeSolver();
        treeSolver.init(spanningTree);
    }

    // solve for x in Lx = b, with default parameters
    @Override
    public double[] solve(double[] b) {
        solve_init(b);
        for (int i = 0; i < 100; ++i) {
            solve_iter();
        }
        return solve_return();
    }

    // find feasible flow on LSST
    // after calling this once, you can call solve_iteration() many times
    public void solve_init(double[] b) {
        currentB = b;
        currentFlow = treeSolver.solveFlow(b);
    }

    // improve the current flow on a cycle induced by an off-tree edge
    public void solve_iter() {
        int e = edgeSampler.next();
    }

    // return the answer, given the flow state
    public double[] solve_return() {
        double[] voltages = new double[spanningTree.nv];

        // build voltage vector from bottom up
        for (int i = 1; i < spanningTree.nv; i++) {
            int v = order[i];
            int parent = spanningTree.getNode(v).getParent().getId();
            double len = spanningTree.getNode(parent).getLength();

            // V = IR
            voltages[v] = voltages[parent] - currentFlow[v]/len;
        }

        // subtract mean voltage
        double meanVoltage = 0;
        for (int i = 0; i < spanningTree.nv; i++) {
            meanVoltage += voltages[i];
        }
        meanVoltage /= spanningTree.nv;
        for (int i = 0; i < spanningTree.nv; i++) {
            voltages[i] -= meanVoltage;
        }

        return voltages;
    }
}
