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
import lapsolver.algorithms.HighStretchSampler;
import lapsolver.algorithms.LDLDecomposition;
import lapsolver.algorithms.Stretch;
import lapsolver.lsst.SpanningTreeStrategy;
import lapsolver.lsst.StarDecompositionTree;
import lapsolver.util.GraphUtils;
import lapsolver.util.LinearAlgebraUtils;
import lapsolver.util.TreeUtils;

import java.util.ArrayList;
import java.util.Random;

import static lapsolver.algorithms.GraphVertexRemoval.AnswerPair;
import static lapsolver.algorithms.LDLDecomposition.ReturnPair;
import static lapsolver.algorithms.Stretch.StretchResult;

public class KMPSolver extends Solver {
    public SpanningTreeStrategy treeStrategy;
    public Solver baseCaseSolver;

    public KMPSolver childSolver;
    public ConjugateGradientSolver recursiveSolver;
    public double tolerance;
    public int maxIters;
    public boolean watch;

    public Graph reducedGraph;
    public double[] reducedD, ldlDiag;
    public Tree spanningTree;
    public Graph sparsifier;

    public AnswerPair gvrPair;
    public ReturnPair ldlPair;
    public int[] gvrInversePerm;

    // Initialize solver with a spanning tree strategy and a solver to run at the bottom level
    public KMPSolver(SpanningTreeStrategy treeStrategy, Solver baseCaseSolver, int maxIters, double tolerance, boolean watch) {
        this.treeStrategy = treeStrategy;
        this.baseCaseSolver = baseCaseSolver;
        this.maxIters = maxIters;
        this.tolerance = tolerance;
        this.watch = watch;
    }

    // Use PCG and StarDecompositionTree strategies
    public KMPSolver(int maxIters, double tolerance, boolean watch) {
        this(new StarDecompositionTree(), new ConjugateGradientSolver(100, 1e-14), maxIters, tolerance, watch);
    }

    public KMPSolver(SpanningTreeStrategy strat) {
        this(strat, new ConjugateGradientSolver(100, 1e-14), 1000, 1e-8, false);
    }

    // vanilla default parameters
    public KMPSolver() {
        this(new StarDecompositionTree(), new ConjugateGradientSolver(100, 1e-14), 1000, 1e-8, false);
    }

    public void init(Graph graph, double[] d, int maxLevels) {
        System.out.println("INIT: n=" + graph.nv + ", m=" + graph.ne);

        this.graph = graph;
        this.d = d;

        eliminate();

        if (maxLevels == 1 || reducedGraph.nv < 500) {
            childSolver = null;
            baseCaseSolver.init(reducedGraph, reducedD);
        } else {
            GraphUtils.reciprocateWeights(reducedGraph);
            spanningTree = treeStrategy.getTree(reducedGraph);
            GraphUtils.reciprocateWeights(reducedGraph);

            sparsifier = sparsify(reducedGraph, spanningTree);

            childSolver = new KMPSolver(treeStrategy, baseCaseSolver, 5, 0, false);
            childSolver.init(sparsifier, reducedD, maxLevels - 1);

            recursiveSolver = new ConjugateGradientSolver(childSolver, maxIters, tolerance, watch);
            recursiveSolver.init(reducedGraph, reducedD);
        }
    }

    public void init(Graph graph, double[] d) {
        init(graph, d, -1);
    }

    // eliminate low-degree vertices
    public void eliminate() {
        // prepare graph for elimination
        GraphVertexRemoval gvr = new GraphVertexRemoval(graph);
        gvrPair = gvr.solve();
        gvrInversePerm = new int[graph.nv];
        for (int i = 0; i < graph.nv; i++) {
            gvrInversePerm[gvrPair.permutation[i]] = i;
        }

        Graph permutedGraph = GraphUtils.permuteGraph(graph, gvrPair.permutation);
        double[] permutedDiag = LinearAlgebraUtils.applyPerm(gvrPair.permutation, d);

        // perform elimination
        LDLDecomposition ldl = new LDLDecomposition(permutedGraph, permutedDiag);
        ldlPair = ldl.solve(gvrPair.numRemoved);

        // get diagonal from LDL
        ldlDiag = new double[graph.nv];
        for (int i = 0; i < ldlPair.D.ne; i++) {
            if (ldlPair.D.u[i] == ldlPair.D.v[i]) {
                ldlDiag[ldlPair.D.u[i]] += ldlPair.D.weight[i];
            }
        }

        // get new graph + diag from LDL
        reducedGraph = LDLDecomposition.getReducedGraph(ldlPair.D, gvrPair.numRemoved);
        reducedD = new double[reducedGraph.nv];

        // compute excess diagonals
        System.arraycopy(ldlDiag, gvrPair.numRemoved, reducedD, 0, graph.nv - gvrPair.numRemoved);
        for (int u = 0; u < reducedGraph.nv; u++) {
            for (int i = 0; i < reducedGraph.deg[u]; i++) {
                reducedD[u] -= reducedGraph.weights[u][i];
            }
        }
    }

    public double[] solve(double[] b) {
//        System.out.println("SOLVE: n=" + graph.nv + ", m=" + graph.ne);

        double[] outerB = ldlPair.L.applyLInv(LinearAlgebraUtils.applyPerm(gvrPair.permutation, b));
        double[] innerB = new double[reducedGraph.nv];
        System.arraycopy(outerB, gvrPair.numRemoved, innerB, 0, innerB.length);

        double[] innerX;

        if (childSolver == null) {
            // we are at the bottom level
            innerX = baseCaseSolver.solve(innerB);
        }
        else {
            innerX = recursiveSolver.solve(innerB);
        }

        double[] outerX = new double[graph.nv];

        for (int i = 0; i < gvrPair.numRemoved; i++) {
            outerX[i] = outerB[i] / ldlDiag[i];
        }

        System.arraycopy(innerX, 0, outerX, gvrPair.numRemoved, graph.nv - gvrPair.numRemoved);

        return LinearAlgebraUtils.applyPerm(gvrInversePerm, ldlPair.L.applyLTransInv(outerX));
    }

    //Construct a preconditioner for graph
    public static Graph sparsify(Graph graph, Tree spanningTree) {
        EdgeList offEdges;

        //Get off-tree edges, find stretches
        offEdges = TreeUtils.getOffTreeEdges(graph, spanningTree);

        GraphUtils.reciprocateWeights(graph);
        StretchResult stretch = Stretch.compute(graph, spanningTree, offEdges);
        GraphUtils.reciprocateWeights(graph);

        //Blow up graph by 4 * avgstretch * log(numRemoved)
        double k = 4. * (stretch.total / (offEdges.ne + 1) * (Math.log(graph.nv) + 1)) *
                (stretch.total / (offEdges.ne + 1) * (Math.log(graph.nv) + 1)) + 1;

        k = 1;
        Graph blownUpGraph = blowUpTreeEdges(graph, spanningTree, k);

        // find stretches in blown-up graph
        GraphUtils.reciprocateWeights(blownUpGraph);
        StretchResult blownUpStretch = Stretch.compute(blownUpGraph, spanningTree, offEdges);
        GraphUtils.reciprocateWeights(blownUpGraph);

        //Expect to grab q = O(m / log(m)) edges
        double q = 10. * graph.ne / Math.pow( Math.log(graph.ne), 2 );

        //Assign p_e = stretch(e) / (total stretch)
        double[] p = blownUpStretch.allStretches.clone();
        for (int i = 0; i < offEdges.ne; i++) {
            p[i] = q * p[i] / blownUpStretch.total;
            if (p[i] > 1) p[i] = 1;
        }

        //nope. sample n/4 edges deterministically
         q = graph.nv / 4;
         p = HighStretchSampler.compute(blownUpStretch.allStretches, (int)q);

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
            sparsifierEdges.weight[index] = 1 / (spanningTree.weight[u] * k);
            index++;
        }

        //Add off-tree edges
        for (int i : edgesToAdd) {
            sparsifierEdges.u[index] = offEdges.u[i];
            sparsifierEdges.v[index] = offEdges.v[i];
            sparsifierEdges.weight[index] = offEdges.weight[i];// / p[i];
            index++;
        }

        System.out.println("Stretch: " + stretch.total + " -> " + blownUpStretch.total);
        System.out.println("E[q] = " + q + ", q = " + edgesToAdd.size());

        // checkSparsifier(graph, new Graph(sparsifierEdges));

        //Build sparsified graph
        return new Graph(sparsifierEdges);
    }

    // given graph and tree, blow up edge weights (not lengths!!) of tree edges in G
    public static Graph blowUpTreeEdges(Graph graph, Tree spanningTree, double k) {
        Graph auxGraph = new Graph(graph);

        // blow up in graph
        for (int u = 0; u < auxGraph.nv; u++) {
            for (int i = 0; i < auxGraph.deg[u]; i++) {
                int v = auxGraph.nbrs[u][i];

                if (spanningTree.parent[u] == v || spanningTree.parent[v] == u) {
                    auxGraph.weights[u][i] *= k;
                }
            }
        }

        // blow up in tree
        for (int i = 0; i < spanningTree.nv; i++) {
            spanningTree.weight[i] /= k;
        }

        return auxGraph;
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
