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
import lapsolver.util.TreeUtils;

import java.util.*;

import static lapsolver.algorithms.GraphVertexRemoval.AnswerPair;
import static lapsolver.algorithms.LDLDecomposition.ReturnPair;
import static lapsolver.algorithms.Stretch.StretchResult;

public class KMP2Solver extends Solver {
    private static final double Cs = 4.0;
    private static final int cStop = 500;
    private static final double kappaC = 1000;

    private final SpanningTreeStrategy treeStrategy;
    private final double failureProbability;
    private static Solver baseCaseSolver = new ConjugateGradientSolver(1000, 1e-9);
    private double[] del;

    public List<ChainEntry> chain;

    public KMP2Solver(SpanningTreeStrategy strategy, double p) {
        treeStrategy = strategy;
        failureProbability = p;
    }

    @Override
    public void init(Graph graph, double[] delta) {
        del = delta;
        chain = buildChain(graph, failureProbability, del);
    }

    @Override
    public double[] solve(double[] b) {
        if (b.length < cStop) {
            baseCaseSolver.init(chain.get(0).graph, del);
            return baseCaseSolver.solve(b);
        }
        return new double[0];
    }

    /**
     * This is the SAMPLE procedure from KMP2.
     *
     * @param edges The edges to sample
     * @param pp    The frequencies of each edge in edges. Indices must correspond to those in edges.
     * @param xi    The error with which to sample
     * @return the chosen edges
     */
    public EdgeList sample(EdgeList edges, StretchResult pp, double xi) {
        // Step 1:  t := \sum_e p'_e
        double t = pp.total;

        // Step 2:  q := C_s t \log t \log (1/\epsilon)
//        final double q = Cs * t * Math.log(t) * Math.log(1 / xi);
        final double q = Cs * t / Math.log(t) / Math.log(t);

        // Step 3:  p_e := p'_e / t
        double[] p = new double[pp.allStretches.length];
        for (int i = 0; i < p.length; i++)
            p[i] = pp.allStretches[i] / t;

        ArrayList<Integer> edgesToAdd = new ArrayList<>((int) q);

        // Using Dan's sampling method.
        for (int i = 0; i < p.length; i++)
            if (Math.random() < p[i])
                edgesToAdd.add(i);

        EdgeList sampledEdges = new EdgeList(edgesToAdd.size());
        for (int i = 0; i < edgesToAdd.size(); i++) {
            int e = edgesToAdd.get(i);
            // Step 7: Add sample of e, l to L_e with weight w'_l  = w_e/(p_e q)
            sampledEdges.u[i] = edges.u[e];
            sampledEdges.v[i] = edges.v[e];
            sampledEdges.weight[i] = edges.weight[e] / (p[e] * q);
        }

        return sampledEdges;
    }

    /**
     * This is IncrementalSparsify from KMP2
     *
     * @param inG the graph to sparsify (resistances)
     * @param inT a low-stretch spanning tree for the graph (resistances)
     * @param kap kappa, the condition number
     * @param xi  epsilon, the maximum allowable error
     * @return the sparsified graph
     */
    public Graph incrementalSparsify(Graph inG, Tree inT, double kap, double xi) {
        if (kap <= 1.0)
            throw new IllegalArgumentException("kappa (" + kap + ") must be > 1");
        if (xi <= 0.0 || xi >= 1.0)
            throw new IllegalArgumentException("0 < xi < 1");
        if (inG.nv != inT.nv)
            throw new IllegalArgumentException("Graph (" + inG.nv + ") and Tree (" + inT.nv
                    + ") must have same number of vertices!");

        // Step 1: Compute stretch_T(G)
        double totalStretch = Stretch.compute(inG, inT, TreeUtils.getOffTreeEdges(inG, inT)).total;

        // Copy the original graph, so we can blow it up
        Tree tPrime = new Tree(inT);

        // Step 2: if |stretch_T(G)| <= 1
        if (totalStretch <= 1) {
            // Step 3: return 2T (resistances -> divide)
            for (int i = 0; i < tPrime.weight.length; i++)
                tPrime.weight[i] /= 2;
            return new Graph(tPrime);
        } // Step 4: end if

        // Step 5: T' := kT - divide because they're resistances
        for (int i = 0; i < tPrime.weight.length; i++)
            tPrime.weight[i] /= kap;

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

        // Step 7: \hat{t} := |stretch_{T'}(G')| = |stretch_{T}(G)| / k
        double tHat = totalStretch / kap;

        // Step 8: t := \hat{t} + n - 1
        final double t = tHat + inT.nv - 1;

        // Implementation step: need off-tree edges:
        final EdgeList offTreeEdges = TreeUtils.getOffTreeEdges(gPrime, tPrime);

        // Step 9: H~ = (V, L~) := SAMPLE(G', stretch_T'(E'), \xi)
        StretchResult stretchResult = Stretch.compute(gPrime, tPrime, offTreeEdges);
        EdgeList hSquiggle = sample(offTreeEdges, stretchResult, xi);

        // Steps 10-12: Check for failure
        final double upperBound = 2 * (tHat / t) * Cs * Math.log(t) * Math.log(1 / xi);
        final double offTreeStretch = stretchResult.total;
        if (offTreeStretch >= upperBound) {
            System.err.println("Error! " + offTreeStretch + " > " + upperBound);
            return null;
        }

        // Step 13: L = L~ - \bigcup_{e \in E_T} L~_{e}
        // since we have no multi edges, this is just the off-tree edges

        // Step 14+15: H := 4(L + 3T') = 4L + 12T'
        int nOffTree = hSquiggle.ne;
        EdgeList tEdges = new EdgeList(tPrime);
        EdgeList H = new EdgeList(nOffTree + tEdges.ne);

        // 4L
        for (int i = 0; i < nOffTree; i++) {
            H.u[i] = hSquiggle.u[i];
            H.v[i] = hSquiggle.v[i];
            H.weight[i] = hSquiggle.weight[i] / 4;
        }

        // 12T'
        for (int i = 0; i < tEdges.ne; i++) {
            H.u[nOffTree + i] = tEdges.u[i];
            H.v[nOffTree + i] = tEdges.v[i];
            H.weight[nOffTree + i] = tEdges.weight[i] / 12;
        }

        return new Graph(H);
    }

    /**
     * BuildChain from KMP2 - Algorithm BuildChain generates the chain of graphs.
     *
     * @param graph The graph to precondition
     * @param p     Failure probability (lower = more success, takes longer)
     * @param delta The
     * @return Chain of Graphs C = {G1, H1, G2, H2, ..., Gd} (no final H!!!)
     */
    public List<ChainEntry> buildChain(Graph graph, double p, double[] delta) {
        // C := 0
        LinkedList<ChainEntry> chain = new LinkedList<>();

        // G1 = G
        Graph g1 = new Graph(graph);

        // T = LowStretchTree(G)
        Tree t = treeStrategy.getTree(graph);

        // H1 = G1 + O~(log^2 n)T
        Graph h1 = new Graph(g1);
        double logSquaredFactor = Math.pow(Math.log(t.nv) * Math.log(Math.log(t.nv)), 2.0);
        int[] parent = t.parent;
        for (int i = 0; i < parent.length; i++)
            if (i != t.parent[i])
                for (int j = 0; j < h1.nbrs[i].length; j++)
                    if (h1.nbrs[i][j] == t.parent[i])
                        h1.weights[i][j] = logSquaredFactor * t.weight[i];

        // G2 = H1
        Graph g2 = new Graph(h1);

        // xi := 2 log n
        double xi = 2 * Math.log(t.nv);

        // Initialize the chain
        chain.add(new ChainEntry(g1, t, logSquaredFactor)); // G
        chain.add(new ChainEntry(h1, t, logSquaredFactor)); // H
        chain.add(new ChainEntry(g2, t, kappaC));           // G

        double[][] deltaRef = new double[1][];
        deltaRef[0] = delta;

        ChainEntry chainEnd = chain.getLast();
        while (chainEnd.graph.nv > cStop) {
            Graph hGraph = incrementalSparsify(chainEnd.graph, chainEnd.tree, chainEnd.kappa, p * xi);
            if (hGraph == null) break; // Failed to sparsify further
            chain.add(new ChainEntry(hGraph, chainEnd.tree, chainEnd.kappa)); // H
            chain.add(greedyElimination(hGraph, chainEnd.tree, deltaRef)); // G
            chainEnd = chain.getLast();
        }

        return chain;
    }

    public ChainEntry greedyElimination(Graph graph, Tree tree, double[][] deltaRef) {
        AnswerPair gvr = new GraphVertexRemoval(graph).solve();
        Graph permutedGraph = GraphUtils.permuteGraph(graph, gvr.permutation);

        deltaRef[0] = applyPerm(gvr.permutation, deltaRef[0]);
        ReturnPair ldl = new LDLDecomposition(permutedGraph, deltaRef[0]).solve(gvr.numRemoved);

        Tree permutedTree = TreeUtils.permuteTree(tree, gvr.permutation);
        Tree updatedTree = updateTree(permutedGraph, permutedTree, ldl.L, gvr.numRemoved);

        Graph reducedGraph = LDLDecomposition.getReducedGraph(permutedGraph, ldl.D, gvr.numRemoved);

        deltaRef[0] = Arrays.copyOfRange(deltaRef[0], gvr.numRemoved, gvr.numRemoved + reducedGraph.nv);

        ChainEntry result = new ChainEntry(reducedGraph, updatedTree, kappaC);
        result.perm = gvr.permutation;
        return result;
    }

    private static double[] applyPerm(int[] perm, double[] x) {
        double[] xPerm = new double[perm.length];

        for (int i = 0; i < x.length; i++)
            xPerm[i] = x[perm[i]];

        return xPerm;
    }

    /**
     * Removes numRemoved vertices from the Tree as per the procedure in KMP2.
     *
     * @param graph      The graph of which tree is a spanning tree
     * @param tree       The spanning tree to update
     * @param lMatrix    The L-matrix of the graph after LDLDecompostion
     * @param numRemoved The number of degree 1 and 2 vertices removed
     * @return The spanning tree with appropriately-weighted edges
     */
    private Tree updateTree(Graph graph, Tree tree, EdgeList lMatrix, int numRemoved) {
        int[] deg = new int[numRemoved];
        for (int i = 0; i < lMatrix.ne; i++)
            if (lMatrix.u[i] != lMatrix.v[i])
                deg[lMatrix.v[i]]++;

        int index = 0;
        EdgeList updatedTree = new EdgeList(tree.nv - numRemoved);
        for (int i = numRemoved; i < tree.nv; i++)
            // Don't include self-loops
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
                    updatedTree.weight[index++] = (graph.weights[u][0] * graph.weights[u][1])
                            / (graph.weights[u][0] + graph.weights[u][1]);
                }
            }

        // Confirm we added back all the edges
        assert index == (tree.nv - numRemoved - 1);
        return new Tree(updatedTree);
    }

    public static class ChainEntry {
        public Graph graph;
        public Tree tree;
        public double kappa;

        public int[] perm;

        public ChainEntry(Graph graph, Tree tree, double kappa) {
            this.graph = graph;
            this.tree = tree;
            this.kappa = kappa;
            this.perm = new int[graph.nv];
            for (int i = 0; i < perm.length; i++) perm[i] = i;
        }
    }
}
