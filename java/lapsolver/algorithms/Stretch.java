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

public class Stretch {
    private final Graph graph;
    private final Tree spanningTree;

    // constructor: just initialize pointers
    public Stretch (Graph graph, Tree spanningTree) {
        this.graph = graph;
        this.spanningTree = spanningTree;
    }

    // do stretch computation for every edge in the tree
    public double totalStretch() {
        double[][] ijv = graph.toIJV();
        double[] pathlen;
        double total = 0;

        // get tree path lengths
        TreePath tp = new TreePath(spanningTree);
        pathlen = tp.query(ijv[0], ijv[1]);

        // divide by edge length for stretch; accumulate
        for (int i = 0; i < graph.ne; i++) {
            total += pathlen[i] / ijv[2][i];
        }

        return total;
    }

    // get mean stretch of edges
    public double meanStretch() {
        return totalStretch() / graph.ne;
    }
}
