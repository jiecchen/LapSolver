/**
 * @file TreeSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Fri Jun 6 2014
 *
 * A linear-time Laplacian solver for trees.
 * Can take a graph, which it converts into a tree, but also
 * has a method that allows it to take a tree directly.
 */

package lapsolver.solvers;

import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.solvers.Solver;
import lapsolver.util.GraphUtils;
import lapsolver.util.TreeUtils;

public class TreeSolver implements Solver {
    private Tree tree;

    /**
     * initialize solver on a graph (must contain a tree)
     *
     * @param graph the graph to solve
     */
    public void init(Graph graph) {
        tree = GraphUtils.toTree(graph);
    }

    /**
     * use a tree directly (not part of interface)
     *
     * @param tree the tree to solve
     */
    public void init(Tree tree) {
        this.tree = tree;
    }

    /**
     * solve for x in Lx = b by calculating potentials down the tree
     *
     * @param b the right hand side of the equations
     * @return x = L^-1 b
     */
    public double[] solve(double[] b) {
        int[] order = TreeUtils.bfsOrder(tree);

        // compute currents from bottom up
        double[] flowTo = new double[tree.nv]; // total flow to permutation
        double[] flowUp = new double[tree.nv]; // flow along (permutation -> parent)
        for (int i = tree.nv - 1; i >= 1; i--) {
            int v = order[i];
            int parent = tree.parent[v];

            flowUp[v] = flowTo[v] - b[v];
            flowTo[parent] += flowUp[v];
        }

        // compute voltages from top down
        double[] voltages = new double[tree.nv];
        for (int i = 1; i < tree.nv; i++) {
            int v = order[i];
            int parent = tree.parent[v];
            double len = tree.weight[v];

            // V = IR
            voltages[v] = voltages[parent] - flowUp[v] * len;
        }

        // subtract mean voltage
        double meanVoltage = 0;
        for (int i = 0; i < tree.nv; i++) {
            meanVoltage += voltages[i];
        }
        meanVoltage /= tree.nv;
        for (int i = 0; i < tree.nv; i++) {
            voltages[i] -= meanVoltage;
        }

        return voltages;
    }

    /**
     * find the feasible flow that gives i_out = b
     *
     * @param b the right hand side of the equations
     * @return a parent list of potentials
     */
    public double[] solveFlow(double[] b) {
        int[] order = TreeUtils.bfsOrder(tree);

        // compute currents from bottom up
        double[] flowTo = new double[tree.nv]; // total flow to permutation
        double[] flowUp = new double[tree.nv]; // flow along (permutation -> parent)
        for (int i = tree.nv - 1; i >= 1; i--) {
            int v = order[i];
            int parent = tree.parent[v];

            flowUp[v] = flowTo[v] - b[v];
            flowTo[parent] += flowUp[v];
        }

        return flowUp;
    }

}
