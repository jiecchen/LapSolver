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

import java.util.*;

import static lapsolver.algorithms.GraphVertexRemoval.AnswerPair;
import static lapsolver.algorithms.LDLDecomposition.LMatrix;
import static lapsolver.algorithms.LDLDecomposition.ReturnPair;
import static lapsolver.algorithms.Stretch.StretchResult;

public class KMP2Solver extends Solver {
    private static final double Cs = 10.0;
    private static final int cStop = 1000;
    private static final double kappaC = 1.0;
    private static final int minIters = 5;
    private final SpanningTreeStrategy treeStrategy;
    public ConjugateGradientSolver recursiveSolver;
    private int maxIters;
    private double tol;
    public LinkedList<ChainEntry> chain;

    public KMP2Solver(SpanningTreeStrategy strategy) {
        this(strategy, 1000);
    }

    public KMP2Solver(SpanningTreeStrategy strategy, int maxIters) {
        this(strategy, maxIters, 1e-14);
    }

    public KMP2Solver(SpanningTreeStrategy strategy, int maxIters, double tolerance) {
        treeStrategy = strategy;
        this.maxIters = maxIters;
        this.tol = tolerance;
    }

    @Override
    public void init(Graph graph, double[] d) {
        chain = buildChain(this.graph = graph, this.d = d);
        recursiveSolver = buildRecursiveSolver(chain, 0);
    }

    @Override
    public double[] solve(double[] b) {
        return recursiveSolver.solve(b);
    }

    private ConjugateGradientSolver buildRecursiveSolver(List<ChainEntry> chain, int level) {
        final ChainEntry current = chain.get(level);
        final int nIters = (level == 0) ? maxIters : minIters;

        if (level == chain.size() - 1)
            return new ConjugateGradientSolver(nIters, 1e-14).initialize(current.graph, current.delta);

        final ChainEntry next = chain.get(level + 1);
        final int numRemoved = current.graph.nv - next.graph.nv;
        final Solver nextSolver = buildRecursiveSolver(chain, level + 1);

        return new ConjugateGradientSolver(new Solver() {
            private int[] invPerm;

            @Override public void init(Graph graph, double[] d) {
                this.graph = graph;
                this.d = d;
                invPerm = new int[graph.nv];
                for (int i = 0; i < current.sparsifier.nv; i++) invPerm[next.perm[i]] = i;
            }

            @Override public double[] solve(double[] b) {
                double[] outerB = next.lMatrix.applyLInv(LinearAlgebraUtils.applyPerm(next.perm, b));
                double[] innerB = Arrays.copyOfRange(outerB, numRemoved, outerB.length);
                double[] innerX = nextSolver.solve(innerB);
                double[] outerX = new double[graph.nv];
                for (int i = 0; i < numRemoved; i++)
                    outerX[i] = outerB[i] / next.diag[i];
                System.arraycopy(innerX, 0, outerX, numRemoved, graph.nv - numRemoved);
                return LinearAlgebraUtils.applyPerm(invPerm, next.lMatrix.applyLTransInv(outerX));
            }
        }.initialize(current.sparsifier, current.delta), nIters, tol, true).initialize(current.graph, current.delta);
    }

    /**
     * BuildChain from KMP2 - Algorithm BuildChain generates the chain of graphs.
     *
     * @param graph The graph to precondition (weights on edges)
     * @param delta The diagonal vector to add to make the Laplacian non-singular
     * @return Chain of Graphs C = {G1, H1, G2, H2, ..., Gd}
     */
    private LinkedList<ChainEntry> buildChain(Graph graph, double[] delta) {
        // C := 0
        LinkedList<ChainEntry> chain = new LinkedList<>();

        // G1 = G
        Graph g1 = new Graph(graph);

        // Initialize the chain
        GraphUtils.reciprocateWeights(graph);
        chain.add(new ChainEntry(g1, treeStrategy.getTree(graph), d));
        chain.getLast().lMatrix = new LDLDecomposition(g1, d).solve(0).L;
        GraphUtils.reciprocateWeights(graph);

        for (ChainEntry chainEnd = chain.getLast(); chainEnd.graph.nv > cStop; chainEnd = chain.getLast()) {
            chainEnd.sparsifier = incrementalSparsify(chainEnd.graph, chainEnd.tree);
            chain.add(greedyElimination(chainEnd.sparsifier, chainEnd.tree));
        }

        // Don't precondition at the end
        chain.getLast().sparsifier = chain.getLast().graph;

        return chain;
    }

    /**
     * This is IncrementalSparsify from KMP2
     *
     * @param inG the graph to sparsify (weights)
     * @param inT a low-stretch spanning tree for the graph (weights)
     * @return the sparsified graph
     */
    private Graph incrementalSparsify(Graph inG, Tree inT) {
        // Copy the original graph, so we can blow it up
        Tree tPrime = new Tree(inT);

        // Step 6: G' := G + (k-1)T  ie. replace T with T'
        Graph gPrime = new Graph(inG);
        for (int u = 0; u < tPrime.parent.length; u++) {
            int v = tPrime.parent[u]; // Now we have edge (u,v)
            if (u != v)
                for (int iV = 0; iV < gPrime.nbrs[u].length; iV++)
                    if (gPrime.nbrs[u][iV] == v) { // this is our edge -- update both sides
                        final double trWt = tPrime.weight[u] * kappaC;
                        final int iU = gPrime.backInd[u][iV];
                        gPrime.weights[u][iV] = trWt;
                        gPrime.weights[v][iU] = trWt;
                    }
        }

        // Step 5: T' := kT
        for (int i = 0; i < tPrime.weight.length; i++)
            tPrime.weight[i] /= kappaC;

        // Implementation step: need off-tree edges:
        final EdgeList offTreeEdges = TreeUtils.getOffTreeEdges(gPrime, tPrime);

        // Step 9: H~ = (V, L~) := SAMPLE(G', stretch_T'(E'), \xi)
        GraphUtils.reciprocateWeights(gPrime);
        StretchResult stretch = Stretch.compute(tPrime, TreeUtils.getOffTreeEdges(gPrime, tPrime));
        GraphUtils.reciprocateWeights(gPrime);
        EdgeList hSquiggle = sample(offTreeEdges, stretch);

        // Step 14+15: H := 4(L + 3T') = 4L + 12T' Implemented here as L+T'
        int nOffTree = hSquiggle.ne;
        EdgeList tEdges = new EdgeList(tPrime);
        EdgeList H = new EdgeList(nOffTree + tEdges.ne);

        for (int i = 0; i < nOffTree; i++) {
            H.u[i] = hSquiggle.u[i];
            H.v[i] = hSquiggle.v[i];
            H.weight[i] = hSquiggle.weight[i];
        }

        for (int i = 0; i < tEdges.ne; i++) {
            H.u[nOffTree + i] = tEdges.u[i];
            H.v[nOffTree + i] = tEdges.v[i];
            H.weight[nOffTree + i] = 1 / tEdges.weight[i];
        }

        return new Graph(H);
    }

    private ChainEntry greedyElimination(Graph graph, Tree tree) {
        AnswerPair gvr = new GraphVertexRemoval(graph).solve();
        Graph permutedGraph = GraphUtils.permuteGraph(graph, gvr.permutation);

        d = LinearAlgebraUtils.applyPerm(gvr.permutation, d);
        ReturnPair ldl = new LDLDecomposition(permutedGraph, d).solve(gvr.numRemoved);

        Graph reducedGraph = LDLDecomposition.getReducedGraph(ldl.D, gvr.numRemoved);
        Tree reducedTree = updateTree(reducedGraph, TreeUtils.permuteTree(tree, gvr.permutation), gvr.numRemoved);

        d = Arrays.copyOfRange(d, gvr.numRemoved, d.length);

        ChainEntry result = new ChainEntry(reducedGraph, reducedTree, d.clone());
        result.perm = gvr.permutation;
        result.lMatrix = ldl.L;
        result.diag = new double[graph.nv];
        for (int i = 0; i < ldl.D.ne; i++)
            if (ldl.D.u[i] == ldl.D.v[i])
                result.diag[ldl.D.u[i]] += ldl.D.weight[i];

        System.arraycopy(result.diag, gvr.numRemoved, d, 0, reducedGraph.nv);
        for (int u = 0; u < reducedGraph.nv; u++)
            for (int i = 0; i < reducedGraph.deg[u]; i++)
                d[u] -= reducedGraph.weights[u][i];

        return result;
    }

    /**
     * This is the SAMPLE procedure from KMP2.
     *
     * @param edges The edges to sample
     * @param pp    The frequencies of each edge in edges. Indices must correspond to those in edges.
     * @return the chosen edges
     */
    private EdgeList sample(EdgeList edges, StretchResult pp) {
        final double q = Cs * edges.ne / Math.pow(Math.log(edges.ne), 2.0);

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
            sampledEdges.weight[i] = edges.weight[e]; // supposed to scale by 1/p[i]
        }

        return sampledEdges;
    }

    /**
     * Removes numRemoved vertices from the Tree.
     *
     * @param graph      The graph of which tree is a spanning tree
     * @param tree       The spanning tree to update
     * @param numRemoved The number of degree 1 and 2 vertices removed
     * @return The spanning tree with appropriately-weighted edges
     */
    private Tree updateTree(Graph graph, Tree tree, int numRemoved) {
        int[] bfsOrder = TreeUtils.bfsOrder(tree);
        int[] newParent = tree.parent.clone();

        // Remove all the non-root vertices below numRemoved
        for (int u : bfsOrder)
            if (u < numRemoved && newParent[u] != u)
                for (int v : tree.children[u])
                    newParent[v] = newParent[u];

        // Assign first non-eliminated child of root as root
        int newRoot = -1;
        for (int i = numRemoved; i < tree.nv; i++)
            if (newParent[i] == tree.root) {
                if (newRoot == -1)
                    newRoot = i;
                newParent[i] = newRoot;
            }

        int[] reducedParent = Arrays.copyOfRange(newParent, numRemoved, newParent.length);
        double[] reducedWeight = new double[reducedParent.length];
        for (int i = 0; i < reducedParent.length; i++)
            reducedParent[i] -= numRemoved;

        for (int u = 0; u < graph.nv; u++)
            for (int i = 0; i < graph.deg[u]; i++)
                if (reducedParent[u] == graph.nbrs[u][i])
                    reducedWeight[u] = 1 / graph.weights[u][i];

        return new Tree(reducedParent, reducedWeight);
    }

    public static class ChainEntry {
        public final Graph graph;
        public final Tree tree;

        public Graph sparsifier;
        public int[] perm;
        public LMatrix lMatrix = null;

        public double[] delta = null;
        public double[] diag = null;

        public ChainEntry(Graph graph, Tree tree, double[] delta) {
            this.graph = graph;
            this.tree = tree;
            this.perm = new int[graph.nv];
            this.delta = delta;
            for (int i = 0; i < perm.length; i++) perm[i] = i;
        }
    }

}
