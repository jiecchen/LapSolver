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
import lapsolver.algorithms.TarjanLCA;

public class Stretch {
    private WeightedGraph graph;
    private Tree spanningTree;

    public Stretch (WeightedGraph graph, Tree spanningTree) {
        this.graph = graph;
        this.spanningTree = spanningTree;
    }

    public double totalStretch() {
        return 0.0;
    }

    public double meanStretch() {
        return 0.0;
    }

    public double edgeStretch(int u, int v) {
        return 0.0;
    }
}
