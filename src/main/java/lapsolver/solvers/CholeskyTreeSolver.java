/**
 * @file TreeSolver.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jul 23 2014
 *
 * A linear-time exact Laplacian solver for trees, for positive definite systems.
 */

package lapsolver.solvers;

import lapsolver.Graph;
import lapsolver.algorithms.GraphVertexRemoval;
import lapsolver.algorithms.LDLDecomposition;
import lapsolver.util.GraphUtils;
import lapsolver.util.LinearAlgebraUtils;

public class CholeskyTreeSolver extends Solver {
    public GraphVertexRemoval.AnswerPair gvrPair;
    public LDLDecomposition.ReturnPair ldlPair;
    public double[] ldlDiagonal;
    public int[] gvrInversePerm;

    @Override
    public void init(Graph graph, double[] d) {
        this.graph = graph;
        this.d = d;

        // get elimination permutation, check sanity
        GraphVertexRemoval gvr = new GraphVertexRemoval(graph);
        gvrPair = gvr.solve();
        if (gvrPair.numRemoved != graph.nv-1) {
            throw new IllegalArgumentException("CholeskyTreeSolver only takes a tree");
        }

        gvrInversePerm = new int[graph.nv];
        for (int i = 0; i < graph.nv; i++) {
            gvrInversePerm[gvrPair.permutation[i]] = i;
        }

        // perform full Cholesky decomposition
        Graph permutedGraph = GraphUtils.permuteGraph(graph, gvrPair.permutation);
        double[] permutedDiag = LinearAlgebraUtils.applyPerm(gvrPair.permutation, d);
        LDLDecomposition ldl = new LDLDecomposition(permutedGraph, permutedDiag);
        ldlPair = ldl.solve(gvrPair.numRemoved);

        // extract diagonal for easy access
        ldlDiagonal = new double[graph.nv];
        for (int i = 0; i < ldlPair.D.ne; i++) {
            if (ldlPair.D.u[i] != ldlPair.D.v[i]) {
                throw new IllegalArgumentException("Cholesky factorization failed");
            }
            else {
                ldlDiagonal[ldlPair.D.u[i]] += ldlPair.D.weight[i];
            }
        }
    }

    @Override
    public double[] solve(double[] b) {
        double[] reducedB = ldlPair.L.applyLInv(LinearAlgebraUtils.applyPerm(gvrPair.permutation, b));
        for (int i = 0; i < graph.nv; i++) {
            reducedB[i] /= ldlDiagonal[i];
        }
        return LinearAlgebraUtils.applyPerm(gvrInversePerm, ldlPair.L.applyLTransInv(reducedB));
    }
}
