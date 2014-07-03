/**
 * @file CGSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Tue Jul 1 2014
 *
 * A wrapper for the Apache Commons-provided conjugate gradient solver.
 */

package lapsolver.solvers;

import lapsolver.Graph;
import lapsolver.util.GraphOperators;
import org.apache.commons.math3.linear.*;

public class ConjugateGradientSolver extends Solver {
    public GraphOperators.LaplacianOperator lapOperator;
    public ConjugateGradient conjugateGradient;
    public RealLinearOperator preconditioner;

    public int maxIters;
    public double tolerance;

    public ConjugateGradientSolver(RealLinearOperator preconditioner, int maxIters, double tolerance) {
        this.preconditioner = preconditioner;
        this.maxIters = maxIters;
        this.tolerance = tolerance;
    }

    // don't use a preconditioner
    public ConjugateGradientSolver(int maxIters, double tolerance) {
        this (null, maxIters, tolerance);
    }

    public void init (Graph graph, double[] d) {
        this.graph = graph;
        this.d = d;
        lapOperator = new GraphOperators.LaplacianOperator(graph, d);
        conjugateGradient = new ConjugateGradient(maxIters, tolerance, true);
    }

    public double[] solve (double[] b) {
        RealVector sol = conjugateGradient.solve(lapOperator, preconditioner, new ArrayRealVector(b));
        return sol.toArray();
    }

}
