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

public class TreeSolver implements Solver {
    private Tree tree;

    // initialize solver on a graph (just converts to a tree)
    public void init(Graph graph) {
        tree = GraphUtils.toTree(graph);
    }

    // use a tree directly (not part of interface)
    public void init(Tree tree) {
        this.tree = tree;
    }

    // solve for x in Lx = b
    public double[] solve(double[] b) {

    }
}
