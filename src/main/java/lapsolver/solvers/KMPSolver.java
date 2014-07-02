/**
 * @file KMPSolver.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Fri Jun 20 2014
 *
 * The recursive preconditioning solver from KMP1.
 */

package lapsolver.solvers;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.GraphVertexRemoval;
import lapsolver.algorithms.LDLDecomposition;
import lapsolver.algorithms.Stretch;
import lapsolver.lsst.SpanningTreeStrategy;
import lapsolver.util.GraphUtils;
import lapsolver.util.TreeUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import static lapsolver.algorithms.GraphVertexRemoval.AnswerPair;
import static lapsolver.algorithms.LDLDecomposition.ReturnPair;
import static lapsolver.algorithms.Stretch.StretchResult;

public class KMPSolver {
    private SpanningTreeStrategy treeStrategy;

    //Initialize solver with a spanning tree strategy
    public KMPSolver(SpanningTreeStrategy treeStrategy) {
        this.treeStrategy = treeStrategy;
    }

    public EdgeList LL;
    public EdgeList DD;
    public Graph innerGraph;

    //Initialize solver on a particular graph, and perform preprocessing
    public double[] solve(Graph graph, double b[], int level, double[] addDiag) {
        System.out.println(graph.nv + " " + graph.ne);
        innerGraph = graph;

        if (graph.nv < 500)
            return solveBaseCase(graph, b, addDiag);

        //Build the preconditioner for the given graph
        Graph sparsifier = buildPreconditioner(graph, treeStrategy.getTree(graph));

        //Create the graph that will be used in the recursion (laplacian given by D matrix)
        GraphVertexRemoval gvrElement = new GraphVertexRemoval(sparsifier);
        AnswerPair gvr = gvrElement.solve();

        b = applyPerm(gvr.permutation, b);
        addDiag = applyPerm(gvr.permutation, addDiag);
        Graph permSparsifier = GraphUtils.permuteGraph(sparsifier, gvr.permutation);

        LDLDecomposition ldlElement = new LDLDecomposition(permSparsifier, addDiag);
        ReturnPair ldl = ldlElement.solve(gvr.numRemoved);

        LL = new EdgeList(ldl.L);
        DD = new EdgeList(ldl.D);

        b = LDLDecomposition.applyLInv(ldl.L, b);

        // grab diagonal elements from D matrix
        double[] diagD = new double[graph.nv];
        for (int i = 0; i < ldl.D.ne; i++) {
            if (ldl.D.u[i] == ldl.D.v[i]) {
                diagD[ldl.D.u[i]] += ldl.D.weight[i];
            }
        }

        Graph reducedSparsifier = buildReducedSparsifier(graph, gvr, ldl);

        // for eliminated vertices, x[i] = b[i]/D[i,i]
        double[] x = new double[graph.nv];
        for (int i = 0; i < gvr.numRemoved; i++) {
            x[i] = b[i] / diagD[i];
        }

        double[] smallb = new double[graph.nv - gvr.numRemoved];
        System.arraycopy(b, gvr.numRemoved, smallb, 0, smallb.length);

        double[] smallAddDiag = new double[graph.nv - gvr.numRemoved];
        System.arraycopy(addDiag, gvr.numRemoved, smallAddDiag, 0, smallAddDiag.length);

        double[] KMPx = solve(reducedSparsifier, smallb, level - 1, smallAddDiag);
        System.arraycopy(KMPx, 0, x, gvr.numRemoved, KMPx.length);

        x = LDLDecomposition.applyLTransInv(ldl.L, x);

        int[] inversePerm = new int[graph.nv];
        for (int i = 0; i < graph.nv; i++)
            inversePerm[gvr.permutation[i]] = i;

        return applyPerm(inversePerm, x);
    }

    public static double[] applyPerm(int[] perm, double[] x) {
        double[] answer = new double[perm.length];

        for (int i = 0; i < x.length; i++)
            answer[i] = x[perm[i]];

        return answer;
    }

    //Construct a preconditioner for graph
    public static Graph buildPreconditioner(Graph graph, Tree spanningTree) {
        EdgeList offEdges;

        //Get off-tree edges, find stretches
        offEdges = TreeUtils.getOffTreeEdges(graph, spanningTree);

        GraphUtils.reciprocateWeights(graph);
        StretchResult stretch = Stretch.compute(graph, spanningTree, offEdges);
        GraphUtils.reciprocateWeights(graph);

        //Blow up graph by 4 * avgstretch * log(numRemoved)
        double k = 4. * (stretch.total / (offEdges.ne + 1) * (Math.log(graph.nv) + 1)) *
                        (stretch.total / (offEdges.ne + 1) * (Math.log(graph.nv) + 1)) + 1;
        Graph blownUpGraph = blowUpTreeEdges(graph, spanningTree, k);

        // find stretches in blown-up graph
        GraphUtils.reciprocateWeights(blownUpGraph);
        StretchResult blownUpStretch = Stretch.compute(blownUpGraph, spanningTree, offEdges);
        GraphUtils.reciprocateWeights(blownUpGraph);

        //Expect to grab q = O(m / log(m)) edges
        double q = 10. * graph.ne / Math.log(graph.ne) / Math.log(graph.ne);

        //Assign p_e = stretch(e) / (total stretch)
        double[] p = blownUpStretch.allStretches.clone();
        for (int i = 0; i < offEdges.ne; i++) {
            p[i] = q * p[i] / blownUpStretch.total;
            if (p[i] > 1) p[i] = 1;
        }

        //Sample the edges
        ArrayList<Integer> edgesToAdd = new ArrayList<>();
        for (int i = 0; i < offEdges.ne; i++) {
            if (Math.random() < p[i]) {
                // with probability p[i], take edge i
                edgesToAdd.add(i);
            }
        }

        //Generate EdgeList
        EdgeList sparsifierEdges = new EdgeList(graph.nv - 1 + edgesToAdd.size());
        int index = 0;

        //Add tree edges
        for (int u = 0; u < graph.nv; u++) {
            if (u == spanningTree.root) continue;
            sparsifierEdges.u[index] = u;
            sparsifierEdges.v[index] = spanningTree.parent[u];
            sparsifierEdges.weight[index] = spanningTree.weight[u] * k;
            index++;
        }

        //Add off-tree edges
        for (int i : edgesToAdd) {
            sparsifierEdges.u[index] = offEdges.u[i];
            sparsifierEdges.v[index] = offEdges.v[i];
            sparsifierEdges.weight[index] = offEdges.weight[i] / p[i];
            index++;
        }

        checkSparsifier(graph, new Graph(sparsifierEdges));

        //Build sparsified graph
        return new Graph(sparsifierEdges);
    }

    // given graph and tree, blow up edge weights (not lengths!!) of tree edges in G
    public static Graph blowUpTreeEdges(Graph graph, Tree spanningTree, double k) {
        Graph auxGraph = new Graph(graph);

        for (int u = 0; u < auxGraph.nv; u++)
            for (int i = 0; i < auxGraph.deg[u]; i++) {
                int v = auxGraph.nbrs[u][i];

                if (spanningTree.parent[u] == v || spanningTree.parent[v] == u)
                    auxGraph.weights[u][i] *= k;
            }

        return auxGraph;
    }

    //Construct the graph for the next step of the recursion
    public static Graph buildReducedSparsifier(Graph graph, AnswerPair gvr, ReturnPair ldl) {
        ldl.D = GraphUtils.sanitizeEdgeList(ldl.D);
        ArrayList<Integer> edgesToAdd = new ArrayList<>();

        int index = 0;
        for (int i = 0; i < ldl.D.ne; i++) {
            if (ldl.D.u[i] >= ldl.D.v[i]) continue;
            if (ldl.D.u[i] >= gvr.numRemoved && ldl.D.v[i] >= gvr.numRemoved) {
                edgesToAdd.add(i);
            }
        }

        EdgeList reducedSparsifierEdges = new EdgeList(edgesToAdd.size());

        for (int i : edgesToAdd) {
            reducedSparsifierEdges.u[index] = ldl.D.u[i] - gvr.numRemoved;
            reducedSparsifierEdges.v[index] = ldl.D.v[i] - gvr.numRemoved;
            reducedSparsifierEdges.weight[index] = -ldl.D.weight[i];
            index++;
        }

        return new Graph(reducedSparsifierEdges);
    }

    //Construct the graph laplacian and call the pcg solver
    public double[] solveBaseCase(Graph graph, double[] b, double[] diag) {
//        System.out.println("*********");
//        System.out.println(graph.nv);
//        System.out.println(Arrays.toString(diag));

        ConjugateGradientSolver solver = new ConjugateGradientSolver(1000, 1e-6);
        solver.init(graph, diag);
        return solver.solve(b);
    }

    public static void checkSparsifier(Graph G, Graph H) {
        int count = 0;

        for (int iter = 0; iter < 100; iter++) {
            double[] x = new double[G.nv];
            for (int i = 0; i < x.length; i++) {
                Random r = new Random();
                x[i] = 10 * r.nextDouble();
            }

            double xLgx = 0;
            for (int i = 0; i < G.nv; i++)
                for (int j = 0; j < G.deg[i]; j++) {
                    int u = i;
                    int v = G.nbrs[i][j];
                    double w = G.weights[i][j];

                    xLgx = xLgx + (x[u] - x[v]) * (x[u] - x[v]) * w;
                }

            double xLhx = 0;
            for (int i = 0; i < H.nv; i++)
                for (int j = 0; j < H.deg[i]; j++) {
                    int u = i;
                    int v = H.nbrs[i][j];
                    double w = H.weights[i][j];

                    xLhx = xLhx + (x[u] - x[v]) * (x[u] - x[v]) * w;
                }

            //System.out.println(xLgx + " " + xLhx);
            if (xLhx > xLgx)
                count++;
        }
        System.out.println(count);
    }
}
