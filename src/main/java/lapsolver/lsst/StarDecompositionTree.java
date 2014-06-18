/**
 * @file StarDecompositionTree.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Fri Jun 13 2014
 *
 * An implementation of the star decomposition tree algorithm in EEST05.
 * Just a skeleton for now.
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.ShortestPathTree;

import static lapsolver.lsst.StarDecompositionWorker.Decomposition;

public class StarDecompositionTree implements SpanningTreeStrategy {
    private StarDecompositionWorker starDecompositionWorker;

    @Override
    public Tree getTree(Graph graph) {
        // TODO(Cyril): set beta/n contraction threshold her
        // TODO(Cyril): randomized centre picking strategy?
        starDecompositionWorker = new StarDecompositionWorker(graph);
        return new Tree(getLowStretchTree(graph, 0));
    }

    /**
     * Return the star coloring of a graph, for MATLAB testing.
     *
     * @param graph the input graph
     * @param x0    the source of the inner ball
     * @return an array colored with ball/cone cuts
     */
    public int[] getStarColoring(Graph graph, int x0) {
        starDecompositionWorker = new StarDecompositionWorker(graph);
        starDecompositionWorker.makeStarCut(graph, x0, new ShortestPathTree(graph, x0));
        return starDecompositionWorker.getColors();
    }

    // generate the spanning tree (LowStretchTree in EEST05)
    public EdgeList getLowStretchTree(Graph graph, int x0) {
        // at every level, compute the shortest path tree
        ShortestPathTree sptInstance = new ShortestPathTree(graph, x0);

        // contract small edges
        // TODO(Cyril): implement this! change weights to 0 (don't forget both directions)

        // obtain star coloring
        EdgeList bridges = starDecompositionWorker.makeStarCut(graph, x0, sptInstance);
        if (bridges == null) {
            return new EdgeList(sptInstance.getTree());
        }

        int nColors = bridges.ne + 1;

        // expand contracted edges
        // TODO(Cyril): implement this! change weights back to original

        // TODO(Cyril): uncomment this once we're worthy
        /*
        int[] parent = sptInstance.getParent();
        int[] parentIndex = sptInstance.getParentIndex();

        for (int u : bridges.u) {
            int pos = u;
            while (pos != x0) {
                int ind = parentIndex[pos];

                graph.weights[pos][ind] /= 2;
                graph.weights[parent[pos]][graph.backInd[pos][ind]] /= 2;

                pos = parent[pos];
            }
        }
        */

        // generate induced subgraph
        Decomposition decomp = starDecompositionWorker.splitGraph(graph, nColors);
        EdgeList[] childTreeEdges = new EdgeList[nColors];

        // get trees recursively
        // TODO(Alex): parallelize here
        for (int color = 0; color < nColors; color++) {
            int xi = decomp.newLabels[color == 0 ? x0 : bridges.v[color - 1]];
            childTreeEdges[color] = getLowStretchTree(decomp.subgraphs[color], xi);
        }

        // merge results from child calls
        EdgeList parentTreeEdges = new EdgeList(graph.nv - 1);
        int edgePos = 0;
        for (int color = 0; color < nColors; color++) {
            // merge child with relabeled vertices
            for (int i = 0; i < childTreeEdges[color].ne; i++) {
                parentTreeEdges.u[edgePos] = decomp.originalLabels[color][childTreeEdges[color].u[i]];
                parentTreeEdges.v[edgePos] = decomp.originalLabels[color][childTreeEdges[color].v[i]];
                parentTreeEdges.weight[edgePos] = childTreeEdges[color].weight[i];
                edgePos++;
            }
        }

        // merge bridges into result
        for (int i = 0; i < bridges.ne; i++) {
            parentTreeEdges.u[edgePos] = bridges.u[i];
            parentTreeEdges.v[edgePos] = bridges.v[i];
            parentTreeEdges.weight[edgePos] = bridges.weight[i];
            edgePos++;
        }

        return parentTreeEdges;
    }
}
