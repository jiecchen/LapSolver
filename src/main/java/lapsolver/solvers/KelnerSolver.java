/**
 * @file KelnerSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An implementation of the primal-dual algorithm of Kelner et al.
 */

package lapsolver.solvers;

import java.lang.Math;

import lapsolver.algorithms.DiscreteSampler;
import lapsolver.algorithms.Stretch;
import lapsolver.solvers.Solver;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.EdgeList;
import lapsolver.lsst.SpanningTreeStrategy;
import lapsolver.solvers.kelner.DirectFlowTree;
import lapsolver.solvers.kelner.FlowTree;
import lapsolver.solvers.kelner.KelnerFlowTree;
import lapsolver.util.GraphUtils;
import lapsolver.util.TreeUtils;
import lapsolver.solvers.TreeSolver;

public class KelnerSolver extends Solver {
    public Tree spanningTree;
    private SpanningTreeStrategy treeStrategy;
    private FlowTree flowTree;
    private TreeSolver treeSolver;

    // edge data to be preprocessed
    private EdgeList offEdges;
    private double[] offStretch;
    private double totalStretch;
    private DiscreteSampler edgeSampler;

    // algorithm state
    private int[] order;

    // initialize solver with a spanning tree strategy
    public KelnerSolver(SpanningTreeStrategy treeStrategy) {
        this.treeStrategy = treeStrategy;
    }

    // initialize solver on a particular graph, and perform preprocessing
    @Override
    public void init(Graph graph, double[] d) {
        if (d != null) {
            throw new IllegalArgumentException("Kelner solver only solves singular Laplacian systems");
        }

        this.graph = new Graph(graph);
        GraphUtils.reciprocateWeights(this.graph);

        // compute LSST, cache BFS order
        spanningTree = treeStrategy.getTree(this.graph);
        order = TreeUtils.bfsOrder(spanningTree);

        // get off-tree edges, find stretches, initialize sampler
        offEdges = TreeUtils.getOffTreeEdges(this.graph, spanningTree);
        Stretch.StretchResult stretchResult = Stretch.compute(this.graph, spanningTree, offEdges);
        offStretch = stretchResult.allStretches;
        totalStretch = stretchResult.total + graph.ne;
        edgeSampler = new DiscreteSampler(offStretch);

        // initialize feasible flow finder for spanning tree
        treeSolver = new TreeSolver();
        treeSolver.init(spanningTree);

        // initialize the cycle query data structure
        flowTree = new KelnerFlowTree(spanningTree, offEdges);
    }

    // solve for x in Lx = b, with default parameters
    @Override
    public double[] solve(double[] b) {
        int iters = (int) Math.ceil(totalStretch * Math.log(totalStretch));
        System.out.println("KelnerSolver: " + iters + " iterations");
        return solve (b, iters);
    }

    // solve for x in Lx = b, with number of iterations
    public double[] solve(double[] b, int iters) {
        solve_init(b);
        for (int i = 0; i < iters; ++i) {
            solve_iter();
        }
        return solve_return();
    }

    // find feasible flow on LSST
    // after calling this once, you can call solve_iteration() many times
    public void solve_init(double[] b) {
        flowTree.setTreeFlows(treeSolver.solveFlow(b));
    }

    // improve the current flow on a cycle induced by an off-tree edge
    public void solve_iter() {
        if (offEdges.ne == 0) {
            return;
        }

        int e = edgeSampler.next();
        double drop = flowTree.query(e);
        double resistance = (offStretch[e] + 1) * offEdges.weight[e];

        flowTree.update(e, -drop / resistance);
    }

    // return the answer, given the flow state
    public double[] solve_return() {
        double[] currentFlow = flowTree.getTreeFlows();
        double[] voltages = new double[spanningTree.nv];

        // build voltage vector from bottom up
        for (int i = 1; i < spanningTree.nv; i++) {
            int v = order[i];
            int parent = spanningTree.parent[v];
            double len = spanningTree.weight[v];

            // V = IR
            voltages[v] = voltages[parent] - currentFlow[v]*len;
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

    // get the flow energy (diagnostic)
    public double getEnergy() {
        double energy = 0;

        double[] currentFlow = flowTree.getTreeFlows();
        // sum tree energies
        for (int i = 0; i < currentFlow.length; i++) {
            if (i == spanningTree.root) continue;

            double f = currentFlow[i];
            double r = spanningTree.weight[i];
            energy += f*f*r;
        }

        // sum off-tree energies
        for (int i = 0; i < offEdges.ne; i++) {
            double f = flowTree.offFlow[i];
            double r = flowTree.offEdges.weight[i];
            energy += f*f*r;
        }

        return energy;
    }
}
