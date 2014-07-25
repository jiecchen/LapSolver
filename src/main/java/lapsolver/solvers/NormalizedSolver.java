/**
 * @file NormalizedSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jul 16 2014
 *
 * A wrapper that reduces systems in normalized Laplacians (D^(-1/2) L D^(-1/2))
 * to Laplacian systems.
 *
 * Idea: inv(D^(1/2) L D^(1/2)) = P D^(1/2) inv(L) D^(1/2) P,
 * where P is the projection onto the orthogonal complement of span({ d^(1/2) })
 */

package lapsolver.solvers;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.util.LinearAlgebraUtils;

public class NormalizedSolver extends Solver {
    public Solver solver;
    public double[] sqrtDeg, unitSqrtDeg;

    public NormalizedSolver (Solver solver) {
        this.solver = solver;
    }

    @Override
    public void init(Graph graph, double[] d) {
        if (d != null) {
            throw new IllegalArgumentException("NormalizedSolver only solves singular Laplacian systems");
        }

        // precompute d^(1/2)
        sqrtDeg = new double[graph.nv];
        for (int i = 0; i < graph.nv; i++) {
            sqrtDeg[i] = Math.sqrt( graph.deg[i] );
        }

        unitSqrtDeg = sqrtDeg.clone();
        LinearAlgebraUtils.normalize(unitSqrtDeg);

        // initialize wrapped solver
        this.graph = graph;
        solver.init(graph, null);
    }

    @Override
    public double[] solve(double[] b) {
        double[] b_lap = b.clone();
        applyProjection(b_lap);
        applyDiagonal(b_lap);

        double[] x = solver.solve(b_lap);
        applyDiagonal(x);
        applyProjection(x);

        return x;
    }

    /**
     * Applies projection P in place.
     * @param v The vector to be projected.
     */
    public void applyProjection(double[] v) {
        double dot = LinearAlgebraUtils.dot(v, unitSqrtDeg);

        for (int i = 0; i < graph.nv; i++) {
            v[i] -= dot * unitSqrtDeg[i];
        }
    }

    /**
     * Applies diagonal matrix D^(1/2) in place.
     * @param v The input vector.
     */
    public void applyDiagonal(double[] v) {
        for (int i = 0; i < graph.nv; i++) {
            v[i] *= sqrtDeg[i];
        }
    }
}
