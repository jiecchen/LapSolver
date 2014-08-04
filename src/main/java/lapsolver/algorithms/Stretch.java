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
    /**
     * compute stretches for a list of edges
     * @param spanningTree A spanning tree of the graph (weights are resistances).
     * @param edges The set of edges of which to take the stretch (weights are resistances).
     * @return An array of stretches for each edge, paired with the total stretch.
     */
    public static StretchResult compute(Tree spanningTree, EdgeList edges) {
        // get tree path lengths
        TreePath tp = new TreePath(spanningTree);
        double [] stretches = tp.query(edges.u, edges.v);
        double total = 0;

        // divide by edge length for stretch; accumulate
        for (int i = 0; i < stretches.length; i++) {
            // divide by edge length
            stretches[i] /= edges.weight[i];
            total += stretches[i];
        }

        return new StretchResult(stretches, total);
    }

    // shorthand: compute stretches for all edges
    public static StretchResult compute(Graph graph, Tree spanningTree) {
        return compute(spanningTree, new EdgeList(graph));
    }

    public static class StretchResult {
        StretchResult(double[] allStretches, double total) {
            this.allStretches = allStretches;
            this.total = total;
        }

        public double[] allStretches;
        public double total;
    }
}
