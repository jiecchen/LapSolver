/**
 * @file SolverWrapper.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jul 10 2014
 *
 * A wrapper to perform reductions between positive-definite and Laplacian solvers.
 */

package lapsolver.solvers;

import lapsolver.EdgeList;
import lapsolver.Graph;

public class SolverWrapper extends Solver {
    public Solver solver;
    public boolean isLaplacianSolver;

    public SolverWrapper (Solver solver, boolean isLaplacianSolver) {
        this.solver = solver;
        this.isLaplacianSolver = isLaplacianSolver;
    }

    @Override
    public void init(Graph graph, double[] d) {
        this.graph = graph;
        this.d = d;

        if ((d == null) == isLaplacianSolver) {
            // got what the solver takes; no reduction needed
            solver.init(graph, d);
        }
        else if (d == null) {
            // reduce Laplacian to positive definite

            // copy induced graph on all vertices except last
            int index = 0;
            EdgeList reducedEdges = new EdgeList(graph.ne - graph.deg[graph.nv-1]);
            for (int u = 0; u < graph.nv - 1; u++) {
                for (int i = 0; i < graph.deg[u]; i++) {
                    if (graph.nbrs[u][i] > u) break;
                    reducedEdges.u[index] = u;
                    reducedEdges.v[index] = graph.nbrs[u][i];
                    reducedEdges.weight[index] = graph.weights[u][i];
                    index++;
                }
            }

            // stick the last vertex in the diagonal
            double[] diag = new double[graph.nv - 1];
            for (int i = 0; i < graph.deg[graph.nv - 1]; i++) {
                diag[graph.nbrs[graph.nv - 1][i]] = graph.weights[graph.nv - 1][i];
            }

            solver.init(new Graph(reducedEdges), diag);
        }
        else {
            // reduce positive definite to Laplacian

            // count non-zeroes to get degree of new dummy vertex
            int augmentedDeg = 0;
            for (int i = 0; i < graph.nv; i++) {
                if (d[i] != 0) ++augmentedDeg;
            }

            // copy graph
            EdgeList augmentedEdges = new EdgeList(graph.ne + augmentedDeg);
            EdgeList originalEdges = new EdgeList(graph);

            for (int i = 0; i < graph.ne; i++) {
                augmentedEdges.u[i] = originalEdges.u[i];
                augmentedEdges.v[i] = originalEdges.v[i];
                augmentedEdges.weight[i] = originalEdges.weight[i];
            }

            // add sink node to take diagonal vertices
            int index = graph.ne;
            for (int i = 0; i < graph.nv; i++) {
                if (d[i] != 0) {
                    augmentedEdges.u[index] = graph.nv;
                    augmentedEdges.v[index] = i;
                    augmentedEdges.weight[index] = d[i];
                    index++;
                }
            }

            solver.init(new Graph(augmentedEdges), null);
        }
    }

    @Override
    public double[] solve(double[] b) {
        double[] b_sys;

        if ((d == null) == isLaplacianSolver) {
            // got what the solver takes, just feed b directly
            b_sys = b;
        }
        else if (d == null) {
            // turn into a positive definite system by chopping off last entry
            b_sys = new double[graph.nv - 1];
            for (int i = 0; i < graph.nv - 1; i++) {
                b_sys[i] = b[i];
            }
        }
        else {
            // turn into a Laplacian system by adding an entry
            b_sys = new double[graph.nv + 1];
            double sum = 0;

            // copy entries
            for (int i = 0; i < graph.nv; i++) {
                b_sys[i] = b[i];
                sum += b[i];
            }

            // enforce sum = 0
            b_sys[graph.nv] = -sum;
        }

        double[] x_sys = solver.solve(b_sys);

        if ((d == null) == isLaplacianSolver) {
            // got what the solver takes; return result directly
            return x_sys;
        }
        else if (d == null) {
            // recover Laplacian solution from positive definite
            double[] x = new double[graph.nv];
            double sum = 0;

            // copy entries except last
            for (int i = 0; i < graph.nv-1; i++) {
                x[i] = x_sys[i];
                sum += x[i];
            }

            // set last entry to make sum 0
            x[graph.nv] = -sum;

            return x;
        }
        else {
            // reduce positive definite to Laplacian
            double[] x = new double[graph.nv];
            double mean = 0;

            // copy entries, ignoring last
            for (int i = 0; i < graph.nv; i++) {
                x[i] = x_sys[i];
                mean += x[i];
            }

            // normalize
            mean /= graph.nv;
            for (int i = 0; i < graph.nv; i++) {
                x[i] -= mean;
            }

            return x;
        }
    }
}
