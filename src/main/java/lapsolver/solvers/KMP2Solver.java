/**
 * @author Alex Reinking <alexander.reinking@yale.edu>
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
import lapsolver.util.LinearAlgebraUtils;
import lapsolver.util.TreeUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import static lapsolver.algorithms.GraphVertexRemoval.AnswerPair;
import static lapsolver.algorithms.LDLDecomposition.LMatrix;
import static lapsolver.algorithms.LDLDecomposition.ReturnPair;
import static lapsolver.algorithms.Stretch.StretchResult;

public class KMP2Solver extends Solver {
    private static final double Cs = 10.0;
    private static final int cStop = 500;
    private static final double kappaC = 50.0;
    private static final double tolerance = 1e-14;
    private static final int maxIters = 100;
    private static final int minIters = 5;
    private final SpanningTreeStrategy treeStrategy;
    public List<ChainEntry> chain;
    private int[] iterations;

    public KMP2Solver(SpanningTreeStrategy strategy) {
        treeStrategy = strategy;
    }

    @Override
    public void init(Graph graph, double[] delta) {
        chain = buildChain(this.graph = graph, this.d = delta);
    }

    /**
     * BuildChain from KMP2 - Algorithm BuildChain generates the chain of graphs.
     *
     * @param graph The graph to precondition (weights on edges)
     * @param delta The
     * @return Chain of Graphs C = {G1, H1, G2, H2, ..., Gd}
     */
    private List<ChainEntry> buildChain(Graph graph, double[] delta) {
        // C := 0
        LinkedList<ChainEntry> chain = new LinkedList<>();

        // G1 = G
        Graph g1 = new Graph(graph);

        // T = LowStretchTree(G)
        Tree t = treeStrategy.getTree(graph);

        // G2 = H1 = G1 + O~(log^2 n)T
        Graph h1_g2 = new Graph(g1);
        double logSquaredFactor = Math.pow(Math.log(t.nv) * Math.log(Math.log(t.nv)), 2.0);
        int[] parent = t.parent;
        for (int u = 0; u < parent.length; u++)
            if (u != t.parent[u])
                for (int iV = 0; iV < h1_g2.nbrs[u].length; iV++)
                    if (h1_g2.nbrs[u][iV] == t.parent[u]) {
                        int v = h1_g2.nbrs[u][iV];
                        double adjustedWeight = t.weight[u] * logSquaredFactor;
                        int iU = h1_g2.backInd[u][iV];
                        h1_g2.weights[u][iV] = adjustedWeight;
                        h1_g2.weights[v][iU] = adjustedWeight;
                    }

        // Initialize the chain
        chain.add(new ChainEntry(g1, t, logSquaredFactor, delta));
        chain.getLast().sparsifier = h1_g2;
        chain.add(new ChainEntry(h1_g2, t, logSquaredFactor, delta));
        chain.getLast().lMatrix = new LDLDecomposition(h1_g2, delta).solve(0).L;

        double[][] deltaRef = new double[1][];
        deltaRef[0] = delta;

        ChainEntry chainEnd = chain.getLast();
        while (chainEnd.graph.nv > cStop) {
            Graph hGraph = incrementalSparsify(chainEnd.graph, chainEnd.tree, chainEnd.kappa);
            chainEnd.sparsifier = hGraph;
            chain.add(greedyElimination(hGraph, chainEnd.tree, deltaRef));
            chainEnd = chain.getLast();
        }

        chain.getLast().sparsifier = chain.getLast().graph;

        return chain;
    }

    /**
     * This is IncrementalSparsify from KMP2
     *
     * @param inG the graph to sparsify (weights)
     * @param inT a low-stretch spanning tree for the graph (weights)
     * @param kap kappa, the condition number
     * @return the sparsified graph
     */
    private Graph incrementalSparsify(Graph inG, Tree inT, double kap) {
        if (kap <= 1.0)
            throw new IllegalArgumentException("kappa (" + kap + ") must be > 1");
        if (inG.nv != inT.nv)
            throw new IllegalArgumentException("Graph (" + inG.nv + ") and Tree (" +
                                                       inT.nv + ") must have same number of vertices!");

        // Step 1: Compute stretch_T(G)
        double totalStretch = computeOffTreeStretch(inG, inT).total;

        // Copy the original graph, so we can blow it up
        Tree tPrime = new Tree(inT);

        // Step 2: if |stretch_T(G)| <= 1
        if (totalStretch <= 1) {
            // Step 3: return 2T
            for (int i = 0; i < tPrime.weight.length; i++)
                tPrime.weight[i] *= 2;
            return new Graph(tPrime);
        } // Step 4: end if

        // Step 5: T' := kT
        for (int i = 0; i < tPrime.weight.length; i++)
            tPrime.weight[i] *= kap;

        // Step 6: G' := G + (k-1)T  ie. replace T with T'
        Graph gPrime = new Graph(inG);
        for (int u = 0; u < tPrime.parent.length; u++) {
            int v = tPrime.parent[u]; // Now we have edge (u,v)
            if (u != v)
                for (int iV = 0; iV < gPrime.nbrs[u].length; iV++)
                    if (gPrime.nbrs[u][iV] == v) { // this is our edge -- update both sides
                        final double trWt = tPrime.weight[u];
                        final int iU = gPrime.backInd[u][iV];
                        gPrime.weights[u][iV] = trWt;
                        gPrime.weights[v][iU] = trWt;
                    }
        }

        // Implementation step: need off-tree edges:
        final EdgeList offTreeEdges = TreeUtils.getOffTreeEdges(gPrime, tPrime);

        // Sanity check
        System.out.printf("%f (new) == %f (old) with k = %f\n",
                          computeOffTreeStretch(gPrime, tPrime).total,
                          totalStretch / kap,
                          kap);

        // Step 9: H~ = (V, L~) := SAMPLE(G', stretch_T'(E'), \xi)
        EdgeList hSquiggle = sample(offTreeEdges, computeOffTreeStretch(gPrime, tPrime));

        // Step 14+15: H := 4(L + 3T') = 4L + 12T'
        int nOffTree = hSquiggle.ne;
        EdgeList tEdges = new EdgeList(tPrime);
        EdgeList H = new EdgeList(nOffTree + tEdges.ne);

        // 4L
        for (int i = 0; i < nOffTree; i++) {
            H.u[i] = hSquiggle.u[i];
            H.v[i] = hSquiggle.v[i];
            H.weight[i] = hSquiggle.weight[i] * 4;
        }

        // 12T'
        for (int i = 0; i < tEdges.ne; i++) {
            H.u[nOffTree + i] = tEdges.u[i];
            H.v[nOffTree + i] = tEdges.v[i];
            H.weight[nOffTree + i] = tEdges.weight[i] * 12;
        }

        return new Graph(H);
    }

    private ChainEntry greedyElimination(Graph graph, Tree tree, double[][] deltaRef) {
        AnswerPair gvr = new GraphVertexRemoval(graph).solve();
        Graph permutedGraph = GraphUtils.permuteGraph(graph, gvr.permutation);

        deltaRef[0] = LinearAlgebraUtils.applyPerm(gvr.permutation, deltaRef[0]);
        ReturnPair ldl = new LDLDecomposition(permutedGraph, deltaRef[0]).solve(gvr.numRemoved);

        Graph reducedGraph = LDLDecomposition.getReducedGraph(ldl.D, gvr.numRemoved);

//        Tree permutedTree = TreeUtils.permuteTree(tree, gvr.permutation);
        Tree updatedTree = treeStrategy.getTree(reducedGraph);
//        Tree updatedTree = updateTree(permutedGraph, permutedTree, gvr.numRemoved);

        deltaRef[0] = Arrays.copyOfRange(deltaRef[0], gvr.numRemoved, gvr.numRemoved + reducedGraph.nv);

        ChainEntry result = new ChainEntry(reducedGraph, updatedTree, kappaC, deltaRef[0].clone());
        result.perm = gvr.permutation;
        result.lMatrix = ldl.L;
        result.diag = new double[graph.nv];
        for (int i = 0; i < ldl.D.ne; i++)
            if (ldl.D.u[i] == ldl.D.v[i])
                result.diag[ldl.D.u[i]] += ldl.D.weight[i];

        System.arraycopy(result.diag, gvr.numRemoved, deltaRef[0], 0, graph.nv - gvr.numRemoved);
        for (int u = 0; u < reducedGraph.nv; u++)
            for (int i = 0; i < reducedGraph.deg[u]; i++)
                deltaRef[0][u] -= reducedGraph.weights[u][i];

        return result;
    }

    /**
     * Computes the off-tree stretch of the graph inG with respect to inT
     *
     * @param inG the graph containing the tree
     * @param inT the tree whose stretch we want to compute
     * @return the total stretch and individual edge stretches
     */
    private static StretchResult computeOffTreeStretch(Graph inG, Tree inT) {
        GraphUtils.reciprocateWeights(inG);
        for (int i = 0; i < inT.weight.length; i++) inT.weight[i] = 1 / inT.weight[i];
        StretchResult stretch = Stretch.compute(inG, inT, TreeUtils.getOffTreeEdges(inG, inT));
        for (int i = 0; i < inT.weight.length; i++) inT.weight[i] = 1 / inT.weight[i];
        GraphUtils.reciprocateWeights(inG);
        return stretch;
    }

    /**
     * This is the SAMPLE procedure from KMP2.
     *
     * @param edges The edges to sample
     * @param pp    The frequencies of each edge in edges. Indices must correspond to those in edges.
     * @return the chosen edges
     */
    private EdgeList sample(EdgeList edges, StretchResult pp) {
        final double q = Cs * edges.ne / Math.log(edges.ne) / Math.log(edges.ne);

        double[] p = new double[pp.allStretches.length];
        for (int i = 0; i < p.length; i++)
            p[i] = Math.min(1.0, q * pp.allStretches[i] / pp.total);

        ArrayList<Integer> edgesToAdd = new ArrayList<>((int) q);

        // Using Dan's sampling method.
        for (int i = 0; i < p.length; i++)
            if (Math.random() < p[i])
                edgesToAdd.add(i);

        EdgeList sampledEdges = new EdgeList(edgesToAdd.size());
        for (int i = 0; i < edgesToAdd.size(); i++) {
            int e = edgesToAdd.get(i);
            sampledEdges.u[i] = edges.u[e];
            sampledEdges.v[i] = edges.v[e];
            sampledEdges.weight[i] = edges.weight[e] / p[e];
        }

        return sampledEdges;
    }

    /**
     * Removes numRemoved vertices from the Tree as per the procedure in KMP2.
     *
     * @param graph      The graph of which tree is a spanning tree
     * @param tree       The spanning tree to update
     * @param numRemoved The number of degree 1 and 2 vertices removed
     * @return The spanning tree with appropriately-weighted edges
     */
    private Tree updateTree(Graph graph, Tree tree, int numRemoved) {
        int[] deg = new int[numRemoved];
        for (int i = 0; i < numRemoved; i++)
            deg[i] = graph.nbrs[i].length;

        // Remove the degree 1's
        for (int u = 0; u < deg.length; u++)
            if (deg[u] == 1) {
                int v = graph.nbrs[u][0];
                if (v < numRemoved)
                    deg[v]--;
            }

        int index = 0;
        EdgeList updatedTree = new EdgeList(tree.nv - numRemoved);
        for (int i = numRemoved; i < tree.nv; i++)
            if (tree.parent[i] >= numRemoved && tree.parent[i] != i) {
                updatedTree.u[index] = Math.min(i, tree.parent[i]) - numRemoved;
                updatedTree.v[index] = Math.max(i, tree.parent[i]) - numRemoved;
                updatedTree.weight[index++] = tree.weight[index];
            }

        for (int u = 0; u < deg.length; u++)
            if (deg[u] == 2) {
                int v1 = graph.nbrs[u][0];
                int v2 = graph.nbrs[u][1];

                boolean u_v1 = tree.parent[u] == v1 || tree.parent[v1] == u;
                boolean u_v2 = tree.parent[u] == v2 || tree.parent[v2] == u;

                assert u_v1 || u_v2;

                // Make sure both endpoints were preserved
                if (u_v1 && u_v2 && !(v1 < numRemoved || v2 < numRemoved)) {
                    updatedTree.u[index] = Math.min(v1, v2) - numRemoved;
                    updatedTree.v[index] = Math.max(v1, v2) - numRemoved;
                    double w1 = graph.weights[u][0];
                    double w2 = graph.weights[u][1];
                    updatedTree.weight[index++] = (w1 * w2) / (w1 + w2);
                }
            }

        // Confirm we added back all the edges
        System.out.printf("*** Removed %d vertices from G with %d.\n", numRemoved, tree.nv);
        System.out.printf("*** %d (actual) == %d (expected)\n", index, tree.nv - numRemoved - 1);
        return new Tree(updatedTree);
    }

    private double[] recSolve(double[] b, final List<ChainEntry> chain, final int level) {
        final ChainEntry current = chain.get(level);
        int effectiveIters = iterations[level];

        if (level == chain.size() - 1) {
            Solver baseCaseSolver = new ConjugateGradientSolver(effectiveIters, tolerance);
            baseCaseSolver.init(current.graph, current.delta);
            return baseCaseSolver.solve(b);
        }

        final ChainEntry next = chain.get(level + 1);
        final int numRemoved = current.graph.nv - next.graph.nv;

        Solver preconditioner = new Solver() {
            @Override public void init(Graph graph, double[] d) {
                this.graph = graph;
                this.d = d;
            }

            @Override public double[] solve(double[] b) {
//                if(level == 1) { // Bypass GVR
//                    baseCaseSolver.init(graph, d);
//                    return baseCaseSolver.solve(b);
//                }
                double[] outerB = next.lMatrix.applyLInv(LinearAlgebraUtils.applyPerm(next.perm, b));
                double[] innerB = Arrays.copyOfRange(outerB, numRemoved, outerB.length);
                double[] innerX = recSolve(innerB, chain, level + 1);
                double[] outerX = new double[graph.nv];
                for (int i = 0; i < numRemoved; i++)
                    outerX[i] = outerB[i] / next.diag[i];
                System.arraycopy(innerX, 0, outerX, numRemoved, graph.nv - numRemoved);
                int[] invPerm = new int[graph.nv];
                for (int i = 0; i < graph.nv; i++) invPerm[next.perm[i]] = i;
                return LinearAlgebraUtils.applyPerm(invPerm, next.lMatrix.applyLTransInv(outerX));
            }
        };
        preconditioner.init(current.sparsifier, current.delta);
        ConjugateGradientSolver pcg = new ConjugateGradientSolver(preconditioner, effectiveIters, tolerance);
        pcg.init(current.graph, current.delta);
        return pcg.solve(b);
    }

    public static class ChainEntry {
        public final Graph graph;
        public final Tree tree;
        public final double kappa;

        public Graph sparsifier;
        public int[] perm;
        public LMatrix lMatrix = null;

        public double[] delta = null;
        public double[] diag = null;

        public ChainEntry(Graph graph, Tree tree, double kappa, double[] delta) {
            this.graph = graph;
            this.tree = tree;
            this.kappa = kappa;
            this.perm = new int[graph.nv];
            this.delta = delta;
            for (int i = 0; i < perm.length; i++) perm[i] = i;
        }
    }

    @Override
    public double[] solve(double[] b) {
        System.out.println("chain.size() = " + chain.size());
        iterations = new int[chain.size()];
        for (int i = 0; i < iterations.length; i++) {
            double t = -(chain.size() - 1) / Math.log((double) minIters / maxIters);
            iterations[i] = (int) Math.round(maxIters * Math.exp(-i / t));
        }
        return recSolve(b, chain, 0);
    }


}
