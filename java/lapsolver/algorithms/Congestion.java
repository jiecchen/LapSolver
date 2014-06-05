/**
 * @file Congestion.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @date Wed Jun 4 2014
 *
 * For an input graph G and tree T, returns the tree with congestion values on each edge (.getCongestionTree())
 * Complexity: O( M log N )
 */

package lapsolver.algorithms;


import lapsolver.Tree;
import lapsolver.Graph;
import lapsolver.util.TreeUtils;

public class Congestion {

    // Variables needed for the algorithm
    public final Tree spanningTree;
    public final Graph graph;
    private final int vertexCount;
    private final int edgeCount;
    private final double[] vertexWeights;
    public final Tree congestionTree;

    public Congestion (Graph G, Tree T) {
        this.graph = G;
        this.spanningTree = T;

        this.vertexCount = graph.nv;
        this.edgeCount = graph.ne;

        this.vertexWeights = GetVertexWeights();
        this.congestionTree = GetCongestionTree();
    }

    // For a graph and its tree, gets the lca values for every edge of the graph
    private int[] LcaResult() {
        TarjanLCA LcaSolver = new TarjanLCA(spanningTree);

        int index = 0;
        int[] LeftNeighbors = new int[edgeCount];
        int[] RightNeighbors = new int[edgeCount];

        for (int i = 0; i < vertexCount; i++) {
            int NeighborCount = graph.nbrs[i].length;

            for (int j = 0; j < NeighborCount; j++)
                if (i < graph.nbrs[i][j]) { // Check only for ordered edges
                    LeftNeighbors[index] = i;
                    RightNeighbors[index] = graph.nbrs[i][j];
                    index++;
                }
        }

        return LcaSolver.solve(LeftNeighbors, RightNeighbors);
    }

    // Instead of incrementing each vertex of the tree by the value of the path (a,b) in the original graph,
    // will increment vertex a by weight(a,b), increment vertex b by weight(a,b) and decrement vertex lca(a,b)
    // by 2 * weight(a,b). Special cases when lca(a,b) = a or b are treated surprisingly smooth.
    // This will help us have a complexity of O(M) on this segment of the code.

    // Assigns values to vertexWeights
    private double[] GetVertexWeights() {
        int[] lcaEdgeValues = LcaResult();
        double[] vertexWeights = new double[vertexCount];

        int index = 0;

        for (int i = 0; i < vertexCount; i++) {
            for (int j = 0; j < graph.nbrs[i].length; j++) {
                if (i < graph.nbrs[i][j]) { // Do only for ordered pairs
                    // A (u,v) edge with LCA and Cost

                    int v = graph.nbrs[i][j];
                    int LCA = lcaEdgeValues[index];
                    double Cost = graph.weights[i][j];

                    vertexWeights[i] += Cost;
                    vertexWeights[v] += Cost;
                    vertexWeights[LCA] -= 2 * Cost;

                    index++;
                }
            }
        }

        return vertexWeights;
    }

    // Gets the congestionTree
    private Tree GetCongestionTree() {
        Tree answer = spanningTree;

        int[] bfsOrdering = TreeUtils.bfsOrder(spanningTree);   // Get the BFS ordering

        for (int i = vertexCount - 1; i > 0; i--) {

            int child = bfsOrdering[i];
            int parent = spanningTree.getNode(child).getParent().getId();

            // congestion cost
            answer.getNode(child).setLength(vertexWeights[child] / spanningTree.getNode(child).getLength());

            vertexWeights[parent] += vertexWeights[child];
        }

        return answer;
    }
}
