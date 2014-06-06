/**
 * @file Stretch.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jun 4 2014
 *
 * Compute the stretch of a tree, or an individual edge.
 */

package lapsolver.algorithms;

import lapsolver.Tree;
import lapsolver.Graph;
import lapsolver.EdgeList;

public class Stretch {
    private final Graph graph;
    private final Tree spanningTree;

    private final double totalStretch;
    private final double[] stretches;

    // constructor: just initialize pointers
    public Stretch (Graph graph, Tree spanningTree) {
        this.graph = graph;
        this.spanningTree = spanningTree;

        // get tree path lengths
        EdgeList edges = new EdgeList(graph);
        TreePath tp = new TreePath(spanningTree);
        stretches = tp.query(edges.u, edges.v);
        double total = 0;

        // divide by edge length for stretch; accumulate
        for (int i = 0; i < graph.ne; i++) {
            stretches[i] /= edges.weight[i];
            total += stretches[i];
        }

        totalStretch = total;
    }

    // get list of edge stretches (parallel to edge list)
    public double[] getAll() {
        return stretches;
    }

    // get total stretch of tree
    public double getTotal() {
        return totalStretch;
    }

}
