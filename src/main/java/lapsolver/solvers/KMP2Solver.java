/**
 * @author Alex Reinking <alexander.reinking@yale.edu>
 */
package lapsolver.solvers;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.Stretch;
import lapsolver.lsst.SpanningTreeStrategy;
import lapsolver.util.TreeUtils;
import lapsolver.util.Unboxed;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import static lapsolver.algorithms.Stretch.StretchResult;

public class KMP2Solver implements Solver {
    private static final double Cs = 4.0;
    private static final int cStop = 100;
    private static final double kappaC = 1.0e-3;

    private final SpanningTreeStrategy treeStrategy;
    private final double failureProbability;
    public List<ChainEntry> chain;

    public KMP2Solver(SpanningTreeStrategy strategy, double p) {
        treeStrategy = strategy;
        failureProbability = p;
    }

    @Override
    public void init(Graph graph) {
        chain = buildChain(graph, failureProbability);
    }

    @Override
    public double[] solve(double[] b) {
        return new double[0];
    }

    /**
     * This is the SAMPLE procedure from KMP2.
     *
     * @param edges The edges to sample
     * @param pp    The frequencies of each edge in edges. Indices must correspond to those in edges.
     * @param eps   The error with which to sample
     * @return the chosen edges
     */
    public EdgeList sample(EdgeList edges, StretchResult pp, double eps) {
        // Step 1:  t := \sum_e p'_e
        double t = pp.total;

        // Step 2:  q := C_s t \log t \log (1/\epsilon)
        final double q = Cs * t * Math.log(t) * Math.log(1 / eps);
        assert q > 0;

        // Step 3:  p_e := p'_e / t
        double[] p = pp.allStretches.clone();
        for (int i = 0; i < p.length; i++) p[i] /= t;

        ArrayList<Integer> edgesToAdd = new ArrayList<>((int) q);

        // Using Dan's sampling method.
        for (int i = 0; i < p.length; i++)
            if (Math.random() < p[i])
                edgesToAdd.add(i);

        EdgeList newG = new EdgeList(edgesToAdd.size());
        for (int i = 0; i < edgesToAdd.size(); i++) {
            int e = edgesToAdd.get(i);
            // Step 7: Add sample of e, l to L_e with weight w'_l  = w_e/(p_e q)
            newG.u[i] = edges.u[e];
            newG.v[i] = edges.v[e];
            newG.weight[i] = edges.weight[e] / (p[e] * q);
        }

        return newG;
    }

    /**
     * This is IncrementalSparsify from KMP2
     *
     * @param inG the graph to sparsify
     * @param inT a low-stretch spanning tree for the graph
     * @param kap kappa, the condition number
     * @param eps epsilon, the maximum allowable error
     * @return the sparsified graph
     */
    public Graph incrementalSparsify(Graph inG, Tree inT, double kap, double eps) {
        // Step 1: Compute stretch_T(G)
        double totalStretch = Stretch.compute(inG, inT).total;

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
        int[] parent = tPrime.parent;
        for (int u = 0; u < parent.length; u++) {
            int v = parent[u];
            if (u != v)
                for (int e = 0; e < gPrime.nbrs[u].length; e++)
                    if (e == v) // this is our edge
                        gPrime.weights[u][e] = tPrime.weight[u];
        }

        // Step 7: \hat{t} := |stretch_{T'}(G')|
        double tHat = totalStretch / kap;

        // Step 8: t := \hat{t} + n - 1
        double t = tHat + inG.nv - 1;

        // Implementation step: need off-tree edges:
        final EdgeList offTreeEdges = TreeUtils.getOffTreeEdges(gPrime, tPrime);

        // Step 9: H~ = (V, L~) := SAMPLE(G', stretch_T'(E'), \eps)
        StretchResult stretchResult = Stretch.compute(gPrime, tPrime, offTreeEdges);
        EdgeList hSquiggle = sample(offTreeEdges, stretchResult, eps);

        // Steps 10-12: Check for failure
        final double upperBound = 2 * (tHat / t) * Cs * Math.log(t) * Math.log(1 / eps);
        final double offTreeStretch = stretchResult.total;
        if (offTreeStretch >= upperBound) {
            System.err.println("Error! " + offTreeStretch + " > " + upperBound);
            return null;
        }

        // Step 13: L = L~ - \bigcup_{e \in E_T} L~_{e}
        // since we have no multi edges, this is just the off-tree edges

        // Step 14+15: H := 4(L + 3T') = 4L + 12T'
        int nOffTree = hSquiggle.ne;
        EdgeList H = new EdgeList(nOffTree + tPrime.nv - 1);

        // 4L
        for (int i = 0; i < nOffTree; i++) {
            H.u[i] = hSquiggle.u[i];
            H.v[i] = hSquiggle.v[i];
            H.weight[i] = 4 * hSquiggle.weight[i];
        }

        // 12T'
        for (int i = 0; i < tPrime.nv - 1; i++) {
            if (i == tPrime.parent[i]) continue;
            H.u[nOffTree + i] = i;
            H.v[nOffTree + i] = tPrime.parent[i];
            H.weight[nOffTree + i] = 4 * 3 * tPrime.weight[i];
        }

        return new Graph(H);
    }

    /**
     * BuildChain from KMP2 - Algorithm BuildChain generates the chain of graphs.
     *
     * @param graph The graph to precondition
     * @param p     Failure probability (lower = more success, takes longer)
     * @return Chain of Graphs C = {G1, H1, G2, H2, ..., Gd} (no final H!!!)
     */
    public List<ChainEntry> buildChain(Graph graph, double p) {
        // C := 0
        List<ChainEntry> chain = new LinkedList<>();

        // G1 = G
        Graph g1 = new Graph(graph);

        // T = LowStretchTree(G)
        Tree t = treeStrategy.getTree(graph);

        // H1 = G1 + O~(log^2 n)T
        Graph h1 = new Graph(g1);
        double logSquaredFactor = Math.pow(Math.log(t.nv) * Math.log(Math.log(t.nv)), 2.0);
        int[] parent = t.parent;
        for (int i = 0; i < parent.length; i++) {
            if (i == t.parent[i]) continue;
            for (int j = 0; j < h1.nbrs[i].length; j++)
                if (h1.nbrs[i][j] == t.parent[i])
                    h1.weights[i][j] = logSquaredFactor * t.weight[i];
        }

        // G2 = H1
        Graph g2 = new Graph(h1);

        // eps := 2 log n
        double epsilon = 2 * Math.log(t.nv);

        // Initialize the chain
        chain.add(new ChainEntry(g1, h1, t, logSquaredFactor));
        chain.add(new ChainEntry(g2, null, t, kappaC));

        ChainEntry chainEnd = chain.get(chain.size() - 1);
        while (chainEnd.gGraph.nv > cStop) {
            chainEnd.hGraph = incrementalSparsify(chainEnd.gGraph, chainEnd.tree, chainEnd.kappa, p * epsilon);
            chain.add(greedyElimination(chainEnd.hGraph, chainEnd.tree));
            chainEnd = chain.get(chain.size() - 1);
        }

        return chain;
    }

    /*
     * TODO: Figure out what the heck to do with GraphVertexRemoval / LDLDecomposition... they're unintelligible!
     */
    public ChainEntry greedyElimination(Graph graph, Tree tree) {
        // Convert to LIL format
        ArrayList<LinkedList<Neighbor>> lilGraph = new ArrayList<>(graph.nv);
        for (int i = 0; i < graph.nv; i++)
            lilGraph.add(new LinkedList<Neighbor>());
        for (int u = 0; u < lilGraph.size(); u++)
            for (int i = 0; i < graph.nbrs[u].length; i++)
                lilGraph.get(u).add(new Neighbor(graph.nbrs[u][i], graph.weights[u][i]));

        ArrayList<Integer> newParents = new ArrayList<>(tree.parent.length);
        for (int i : tree.parent) newParents.add(i);

        ArrayList<Double> newTreeWeights = new ArrayList<>(tree.weight.length);
        for (double i : tree.weight) newTreeWeights.add(i);

        boolean removedOne;
        do {
            removedOne = false;
            // Remove degree 1
            for (int u = 0; u < lilGraph.size(); u++)
                if (lilGraph.get(u).size() == 1) {
                    for (int i = 0; i < newParents.size(); i++)
                        if (newParents.get(i) > u)
                            newParents.set(i, newParents.get(i) - 1);

                    newParents.remove(u);
                    newTreeWeights.remove(u);
                    removeVertex(lilGraph, u);
                    removedOne = true;
                }

            // Remove degree 2
            for (int curVtx = 0; curVtx < lilGraph.size(); curVtx++) {
                LinkedList<Neighbor> vertex = lilGraph.get(curVtx);
                if (vertex.size() == 2) {
                    Neighbor u1 = vertex.get(0);
                    Neighbor u2 = vertex.get(1);
                    double weight = 1. / (1. / u1.weight + 1. / u2.weight);

                    lilGraph.get(u1.v).add(new Neighbor(u2.v, weight));
                    lilGraph.get(u2.v).add(new Neighbor(u1.v, weight));

                    boolean v_u1_exists = newParents.get(curVtx) == u1.v || newParents.get(u1.v) == curVtx;
                    boolean v_u2_exists = newParents.get(curVtx) == u2.v || newParents.get(u2.v) == curVtx;

                    // If both edges are in the tree
                    if (v_u1_exists && v_u2_exists) {
                        if (newParents.get(curVtx) == curVtx) { // If I'm the root
                            newParents.set(u1.v, u2.v); // Make u2 the new root
                            newTreeWeights.set(u1.v, weight);
                            newParents.set(u2.v, u2.v);
                        } else if (newParents.get(u1.v) == curVtx) {
                            newParents.set(u1.v, u2.v);
                            newTreeWeights.set(u1.v, weight);
                        } else if (newParents.get(u2.v) == curVtx) {
                            newParents.set(u2.v, u1.v);
                            newTreeWeights.set(u2.v, weight);
                        }
                    } else if (v_u1_exists && newParents.get(curVtx) == curVtx)
                        newParents.set(u1.v, u1.v);
                    else if (v_u2_exists && newParents.get(curVtx) == curVtx)
                        newParents.set(u2.v, u2.v);

                    for (int i = 0; i < newParents.size(); i++)
                        if (newParents.get(i) > curVtx)
                            newParents.set(i, newParents.get(i) - 1);

                    removeVertex(lilGraph, curVtx);
                    newParents.remove(curVtx);
                    newTreeWeights.remove(curVtx);
                    removedOne = true;
                }
            }
        } while (removedOne);

        // Rebuild graph
        int ne = 0;
        for (LinkedList<Neighbor> neighbors : lilGraph) ne += neighbors.size();
        ne /= 2;

        int[] us = new int[ne];
        int[] vs = new int[ne];
        double[] ws = new double[ne];
        int index = 0;

        for (int u = 0; u < lilGraph.size(); u++) {
            for (Neighbor neighbor : lilGraph.get(u))
                if (u < neighbor.v) {
                    us[index] = u;
                    vs[index] = neighbor.v;
                    ws[index++] = neighbor.weight;
                }
        }

        Tree newTree = new Tree(Unboxed.intsToArray(newParents), Unboxed.doublesToArray(newTreeWeights));
        return new ChainEntry(new Graph(us, vs, ws), null, newTree, kappaC);
    }

    private void removeVertex(ArrayList<LinkedList<Neighbor>> lilGraph, int vertex) {
        LinkedList<Neighbor> neighbors = lilGraph.get(vertex);
        for (Neighbor neighbor : neighbors)
            for (Iterator<Neighbor> i = lilGraph.get(neighbor.v).iterator(); i.hasNext(); )
                if (i.next().v == vertex)
                    i.remove();
        lilGraph.remove(vertex);
        for (LinkedList<Neighbor> nbrs : lilGraph)
            for (Neighbor nbr : nbrs)
                if (nbr.v > vertex)
                    nbr.v--;
    }

    private static class Neighbor {
        int v;
        double weight;

        private Neighbor(int v, double weight) {
            this.v = v;
            this.weight = weight;
        }
    }

    public static class ChainEntry {
        public Graph gGraph;
        public Graph hGraph;
        public Tree tree;
        public double kappa;

        private ChainEntry(Graph gGraph, Graph hGraph, Tree tree, double kappa) {
            this.gGraph = gGraph;
            this.hGraph = hGraph;
            this.tree = tree;
            this.kappa = kappa;
        }
    }
}
