/**
 * @file GraphOperators.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jul 2 2014
 *
 * Wrappers for Apache Commons linear operators.
 */

package lapsolver.util;

import lapsolver.Graph;
import lapsolver.LinearOperator;

public class GraphOperators {
    public static LinearOperator laplacian(final Graph graphIn, final double[] dIn) {
        return new LinearOperator() {
            Graph graph = new Graph(graphIn);
            double[] d = dIn.clone();
            @Override
            public double[] apply (double[] x) {
                return applyLaplacian(graph, d, x);
            }
        };
    }

    public static LinearOperator laplacian(final Graph graphIn) {
        return new LinearOperator() {
            Graph graph = new Graph(graphIn);
            @Override
            public double[] apply (double[] x) {
                return applyLaplacian(graph, x, null);
            }
        };
    }

    // apply the Laplacian matrix of a graph to a vector
    public static double[] applyLaplacian(Graph graph, double[] d, double[] x) {
        double[] y = new double[graph.nv];
        for (int u = 0; u < graph.nv; u++) {
            for (int i = 0; i < graph.deg[u]; i++) {
                int v = graph.nbrs[u][i];
                y[u] += graph.weights[u][i] * (x[u] - x[v]);
            }
        }

        if (d != null) {
            for (int u = 0; u < graph.nv; u++) {
                y[u] += d[u] * x[u];
            }
        }

        return y;
    }
}
