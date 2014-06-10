/**
 * @file Petal.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Tue Jun 10 2014
 *
 * Naive implementation of petal finding
 *
 */

package lapsolver.algorithms;

import lapsolver.Graph;

public class Petal {
    private int[] charVec;

    public int[] getCharVec() {
        return charVec;
    }

    public Petal (Graph graph, int x0, int target, double radius) {
        ShortestPathTree tree_x0 = new ShortestPathTree(graph, x0);
        double[] dist_x0 = tree_x0.getDist();
        int[] parent_x0 = tree_x0.getParent();

        int pos = target;
        charVec = new int[graph.nv];

        while (dist_x0[target] - dist_x0[pos] <= radius) {
            double coneRadius = (radius - dist_x0[target] + dist_x0[pos]) / 2;

            // add current cone to set
            Cone cone = new Cone(graph, x0, pos);
            double[] coneMetric = cone.getConeDist();
            for (int i = 0; i < graph.nv; i++) {
                if (coneMetric[i] < coneRadius) {
                    charVec[i] = 1;
                }
            }

            // walk back towards x0
            pos = parent_x0[pos];
            if (pos == x0) break;
        }
    }
}
