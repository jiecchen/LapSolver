/**
 * @file Stretch.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jun 4 2014
 *
 * Algorithms that compute stretches of edges and entire trees.
 */

package lapsolver.algorithms;

import lapsolver.Tree;
import lapsolver.WeightedGraph;
import lapsolver.algorithms.TreePath;

public class Stretch {
    private WeightedGraph graph;
    private Tree spanningTree;

    public Stretch (WeightedGraph graph, Tree spanningTree) {
        this.graph = graph;
        this.spanningTree = spanningTree;
    }

    public double totalStretch() {
        double[][] ijv = graph.toIJV();
        double[] pathlen;
        double total = 0;

        TreePath tp = new TreePath(spanningTree);
        pathlen = tp.query(ijv[0], ijv[1]);

        for (int i = 0; i < graph.ne; i++) {
            total += pathlen[i] / ijv[2][i];
        }

        return total;
    }

    public double meanStretch() {
        return totalStretch() / graph.ne;
    }

    public double edgeStretch(int u, int v) {
        return 0.0;
    }
}
