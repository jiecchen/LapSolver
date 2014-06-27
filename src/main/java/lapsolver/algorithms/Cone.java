/**
 * @file Cone.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Tue Jun 10 2014
 *
 * Naive implementation of cone finding algorithm.
 * Mostly for visualization.
 *
 *
 *
 */

package lapsolver.algorithms;

import lapsolver.Graph;
import lapsolver.algorithms.ShortestPathTree;

public class Cone {
    private double[] coneDist;

    public double[] getConeDist() {
        return coneDist;
    }

    /**
     * @param graph The input graph (weights are resistances).
     * @param x0 The center of the cone.
     * @param x The target of the cone.
     */
    public Cone (Graph graph, int x0, int x) {
        ShortestPathTree tree_x0 = new ShortestPathTree(graph, x0);
        ShortestPathTree tree_x = new ShortestPathTree(graph, x);

        double[] dist_x0 = tree_x0.getDist();
        double[] dist_x = tree_x.getDist();

        coneDist = new double[graph.nv];

        for (int i = 0; i < graph.nv; i++) {
            coneDist[i] = dist_x0[x] + dist_x[i] - dist_x0[i];
        }
    }
}
