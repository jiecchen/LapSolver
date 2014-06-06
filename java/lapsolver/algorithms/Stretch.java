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
    public static StretchResult computeStretch(Graph graph, Tree spanningTree) {
        // get tree path lengths
        EdgeList edges = new EdgeList(graph);
        TreePath tp = new TreePath(spanningTree);
        double [] stretches = tp.query(edges.u, edges.v);
        double total = 0;

        // divide by edge length for stretch; accumulate
        for (int i = 0; i < graph.ne; i++) {
            stretches[i] /= edges.weight[i];
            total += stretches[i];
        }

        return new StretchResult(stretches, total);
    }

    static class StretchResult {
        StretchResult(double[] allStretches, double total) {
            this.allStretches = allStretches;
            this.total = total;
        }

        public double[] allStretches;
        public double total;
    }
}
