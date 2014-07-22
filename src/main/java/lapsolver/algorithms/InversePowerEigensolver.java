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
import lapsolver.util.LinearAlgebraUtils;

public class InversePowerEigensolver {
    public Graph graph;
    public Solver inverseOperator;
    public double[] nullspaceBasis;

    /**
     * Initialize an instance of an eigensolver, with an arbitrary nullspace basis vector.
     */
    public InversePowerEigensolver(Graph graph, Solver inverseOperator, double[] nullspaceBasis) {
        this.graph = graph;
        this.inverseOperator = inverseOperator;
        this.inverseOperator.init(graph);
        this.nullspaceBasis = nullspaceBasis;
    }

    /**
     * Initialize an instance of an eigensolver, with a choice between Laplacian (all-ones)
     * and normalized Laplacian (square root of degrees) nullspace basis vectors
     */
    public InversePowerEigensolver(Graph graph, Solver inverseOperator, boolean isNormalizedLaplacian) {
        this(graph, inverseOperator, null);

        // set up nullspace basis vector
        nullspaceBasis = new double[graph.nv];

        if (isNormalizedLaplacian) {
            for (int i = 0; i < graph.nv; i++) {
                nullspaceBasis[i] = Math.sqrt(graph.deg[i]);
            }
        }
        else {
            for (int i = 0; i < graph.nv; i++) {
                nullspaceBasis[i] = 1;
            }
        }

        LinearAlgebraUtils.normalize(nullspaceBasis);
    }

    /**
     * Initialize an instance of a Laplacian eigensolver (default).
     */
    public InversePowerEigensolver(Graph graph, Solver inverseOperator) {
        this(graph, inverseOperator, false);
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
                if (nullspaceBasis != null) {
                    // orthogonalize with respect to nullspace basis

                    LinearAlgebraUtils.normalize(result[i]);
                    double dot = LinearAlgebraUtils.dot(nullspaceBasis, result[i]);
                    for (int k = 0; k < graph.nv; k++) {
                        result[i][k] -= nullspaceBasis[k] * dot;
                    }
                }

                // orthogonalize with respect to other vectors
                for (int i_prev = 0; i_prev < i; i_prev++) {
                    LinearAlgebraUtils.normalize(result[i]);

                    // subtract projection of (this vec) on (earlier vec)
                    double dot = LinearAlgebraUtils.dot(result[i_prev], result[i]);
                    for (int k = 0; k < graph.nv; k++) {
                        result[i][k] -= result[i_prev][k] * dot;
                    }
                }

                // v <- L^-1 v
                double[] x = inverseOperator.solve(result[i]);
                System.arraycopy(x,0,result[i],0,graph.nv);
            }

            LinearAlgebraUtils.normalize(result[i]);
        }

        return result;
    }

}
