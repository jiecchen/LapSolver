/**
 * @file LaplacianEigs.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Mon Jul 14 2014
 *
 * The power method for Laplacian eigenvector computation.
 * Takes a Laplacian inverse operator.
 */

package lapsolver.algorithms;

import lapsolver.Graph;
import lapsolver.solvers.Solver;

public class LaplacianEigs {
    public Graph graph;
    public Solver inverseOperator;

    /**
     * Initialize an instance of an eigensolver.
     */
    public LaplacianEigs (Graph graph, Solver inverseOperator) {
        this.graph = graph;
        this.inverseOperator = inverseOperator;
        this.inverseOperator.init(graph);
    }

    /**
     * Get the first few eigenvectors of the graph's Laplacian.
     * @param nEigs The number of eigenvectors to compute.
     * @param iters The number of iterations to run the power method.
     * @return A 2D array of eigenvectors, nEigs by graph.nv.
     */
    public double[][] getVectors (int nEigs, int iters) {
        double[][] result = new double[nEigs][graph.nv];

        for (int i = 0; i < nEigs; i++) {
            // start with random vector
            for (int j = 0; j < graph.nv; j++) {
                result[i][j] = Math.random();
            }

            // keep applying L^-1
            for (int j = 0; j < iters; j++) {
                // normalize
                normalize(result[i]);

                // orthogonalize with respect to all-ones
                double sum = 0;
                for (int k = 0; k < graph.nv; k++) {
                    sum += result[i][k];
                }
                for (int k = 0; k < graph.nv; k++) {
                    result[i][k] -= sum / graph.nv;
                }

                // orthogonalize with respect to other vectors
                for (int i_prev = 0; i_prev < i; i_prev++) {
                    normalize(result[i]);

                    // compute (earlier vec) dot (this vec)
                    double dot = 0;
                    for (int k = 0; k < graph.nv; k++) {
                        dot += result[i_prev][k]*result[i][k];
                    }

                    // subtract projection
                    for (int k = 0; k < graph.nv; k++) {
                        result[i][k] -= result[i_prev][k] * dot;
                    }
                }

                // v <- L^-1 v
                double[] x = inverseOperator.solve(result[i]);
                System.arraycopy(x,0,result[i],0,graph.nv);
            }

            normalize(result[i]);
        }

        return result;
    }

    public static void normalize (double[] v) {
        double norm = 0;
        for (int i = 0; i < v.length; i++) {
            norm += v[i] * v[i];
        }
        norm = Math.sqrt(norm);
        for (int i = 0; i < v.length; i++) {
            v[i] /= norm;
        }
    }

}
