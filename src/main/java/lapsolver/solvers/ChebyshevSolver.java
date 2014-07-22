/**
 * @file ChebyshevSolver.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @date Tue Jul 22 2014
 *
 * A wrapper for the Chebyshev Solver. Reference at http://en.wikipedia.org/wiki/Chebyshev_iteration.
 */

package lapsolver.solvers;

import lapsolver.Graph;
import lapsolver.LinearOperator;
import lapsolver.util.GraphOperators;
import lapsolver.util.LinearAlgebraUtils;

public class ChebyshevSolver extends Solver {
    public LinearOperator operator;
    public LinearOperator preconditioner;

    public int maxIters;
    public double tolerance;

    public double[] r;
    public double[] z;
    public double[] p;
    public double[] x;

    public double lMin;
    public double lMax;

    /**
     * Standard constructor.
     * @param preconditioner
     * @param maxIters
     * @param tolerance
     */
    public ChebyshevSolver(Solver preconditioner, int maxIters, double tolerance, double lMin, double lMax) {
        this.preconditioner = preconditioner;
        this.maxIters = maxIters;
        this.tolerance = tolerance;
        this.lMin = lMin;
        this.lMax = lMax;
    }

    /**
     * Constructor with no preconditioner.
     * @param maxIters
     * @param tolerance
     */
    public ChebyshevSolver(int maxIters, double tolerance, double lMin, double lMax) {
        this (null, maxIters, tolerance, lMin, lMax);
    }

    @Override
    public void init (Graph graph, double[] d) {
        this.graph = graph;
        this.d = d;

        // allocate scratch space
        r = new double[graph.nv];
        z = new double[graph.nv];
        p = new double[graph.nv];
        x = new double[graph.nv];

        operator = GraphOperators.laplacian(graph, d);
    }

    @Override
    public double[] solve(double b[]) {
        double d = (lMax + lMin) / 2;
        double c = (lMax - lMin) / 2;
        r = LinearAlgebraUtils.subtract(b, operator.apply(x));

        double alpha = 0;
        double beta = 0;

        for (int step = 1; step <= maxIters; step++) {
            if (preconditioner == null) {
                System.arraycopy(r, 0, z, 0, graph.nv);
            }
            else {
                z = preconditioner.apply(r);
            }

            if (step == 1) {
                p = z.clone();
                alpha = 1 / d;
            }
            else {
                beta = (c * alpha / 2) * (c * alpha / 2);
                alpha = 1 / (d - beta / alpha);
                p = LinearAlgebraUtils.add(z, LinearAlgebraUtils.scale(p, beta));
            }

            x = LinearAlgebraUtils.add(x, LinearAlgebraUtils.scale(p, alpha));
            r = LinearAlgebraUtils.subtract(b, operator.apply(x));

            if (LinearAlgebraUtils.norm(r) < tolerance) {
                break;
            }
        }

        return x;
    }
}
