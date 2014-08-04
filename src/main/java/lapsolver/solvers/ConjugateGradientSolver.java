/**
 * @file ConjugateGradientSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Tue Jul 1 2014
 *
 * A wrapper for the Apache Commons-provided conjugate gradient solver.
 */

package lapsolver.solvers;

import lapsolver.Graph;
import lapsolver.LinearOperator;
import lapsolver.util.GraphOperators;
import lapsolver.util.LinearAlgebraUtils;

import java.util.Arrays;

public class ConjugateGradientSolver extends Solver {
    public LinearOperator preconditioner;
    public LinearOperator operator;

    public int maxIters;
    public double tolerance;
    public boolean watch;

    // state variables for PCG iterations
    public double[] r;
    public double[] z;
    public double[] p;
    public double[] x;

    // introspection state
    public double[] watchNorms; // |Lx - b| at each step
    public int watchIters; // how many iterations actually occurred?

    /**
     * Standard constructor.
     *
     * @param preconditioner The preconditioner M^(-1).
     * @param maxIters       The number of iterations at which to terminate, regardless of residual norm.
     * @param tolerance      The termination threshold for residual norm.
     * @param watch          Do we record norms for posterity?
     */
    public ConjugateGradientSolver(Solver preconditioner, int maxIters, double tolerance, boolean watch) {
        this.preconditioner = preconditioner;
        this.maxIters = maxIters;
        this.tolerance = tolerance;
        this.watch = watch;

        if (watch) {
            this.watchNorms = new double[maxIters];
        }
    }

    public ConjugateGradientSolver initialize(Graph g, double[] d) {
        return (ConjugateGradientSolver) super.initialize(g, d);
    }

    /**
     * Constructor with no preconditioner and no introspection.
     */
    public ConjugateGradientSolver(int maxIters, double tolerance) {
        this(null, maxIters, tolerance, false);
    }

    /**
     * Constructor with preconditioner but no introspection.
     */
    public ConjugateGradientSolver(Solver preconditioner, int maxIters, double tolerance) {
        this(preconditioner, maxIters, tolerance, false);
    }

    /**
     * Constructor with no preconditioner.
     */
    public ConjugateGradientSolver(int maxIters, double tolerance, boolean watch) {
        this(null, maxIters, tolerance, watch);
    }

    @Override
    public void init(Graph graph, double[] d) {
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
    public double[] solve(double[] b) {
        // x0 = 0
        Arrays.fill(x, 0);

        // r0 = b - L * x0 (but x0 = 0)
        System.arraycopy(b, 0, r, 0, graph.nv);

        // z0 = P * r0
        if (preconditioner == null) {
            System.arraycopy(r, 0, z, 0, graph.nv);
        } else {
            z = preconditioner.apply(r);
        }

        // p0 = z0
        System.arraycopy(z, 0, p, 0, graph.nv);

        for (int iter = 0; iter < maxIters; iter++) {
            watchIters = iter + 1;

            // alpha = (r dot z) / (p dot Ap)
            double[] ap = operator.apply(p);
            double alpha = LinearAlgebraUtils.dot(r, z) / LinearAlgebraUtils.dot(p, ap);

            // x = x + alpha * p
            for (int i = 0; i < graph.nv; i++) {
                x[i] += alpha * p[i];
            }

            // save (oldz dot oldr)
            double beta_denom = LinearAlgebraUtils.dot(z, r);

            // r = r - alpha * Ap
            for (int i = 0; i < graph.nv; i++) {
                r[i] -= alpha * ap[i];
            }

            // if residual is small enough, return x
            double residualNorm = LinearAlgebraUtils.norm(r);
            if (watch) {
                watchNorms[iter] = residualNorm;
            }
            if (residualNorm < tolerance) {
                return x;
            }

            // z = P * r
            if (preconditioner == null)
                System.arraycopy(r, 0, z, 0, graph.nv);
            else
                z = preconditioner.apply(r);

            // beta = (z dot r) / (oldz dot oldr)
            double beta = LinearAlgebraUtils.dot(z, r) / beta_denom;

            // p = z + beta*p
            for (int i = 0; i < graph.nv; i++) {
                p[i] = z[i] + beta * p[i];
            }
        }
//        System.out.println("PCG: defaulted");
        return x;
    }

}
