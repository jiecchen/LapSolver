/**
 * @file IncrementalSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Fri Jun 27 2014
 *
 * A "meta-wrapper", which turns a bad solver into an arbitrarily good one.
 */

package lapsolver.solvers;

import lapsolver.Graph;
import lapsolver.util.GraphUtils;

public class IncrementalSolver extends Solver {
    public Solver solver; // The underlying solver.
    public int nIterations;

    public double[] residue;
    public double[] currentX;

    public IncrementalSolver(Solver solver, int nIterations) {
        this.solver = solver;
        this.nIterations = nIterations;
    }

    @Override
    public void init(Graph graph, double[] d) {
        this.graph = graph;
        this.d = d;
        solver.init(graph, d);
    }

    @Override
    public double[] solve(double[] b) {
        solve_init(b);
        for (int i = 0; i < nIterations; i++) {
            solve_iter();
        }
        return solve_return();
    }

    // initialize with a residue of b
    public void solve_init(double[] b) {
        residue = b.clone();
        currentX = new double[graph.nv];
    }

    // perform one iteration of the incremental meta-solving process
    public void solve_iter() {
        double[] dx = solver.solve(residue);

        // update current guess
        for (int i = 0; i < graph.nv; i++) {
            currentX[i] += dx[i];
        }

        // update residue
        double[] ldx = GraphUtils.applyLaPlacian(graph, dx);
        if (d != null) {
            for (int i = 0; i < graph.nv; i++) {
                ldx[i] += d[i] * dx[i];
            }
        }

        for (int i = 0; i < graph.nv; i++) {
            residue[i] -= ldx[i];
        }
    }

    // return the current guess of x
    public double[] solve_return() {
        return currentX;
    }
}
