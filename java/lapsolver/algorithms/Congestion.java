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
    private Tree spanningTree;
    private Graph graph;
    private int vertexCount;
    private int edgeCount;
    private double[] vertexWeights;
    public Tree congestionTree;

    public Congestion (Graph G, Tree T) {
        this.graph = G;
        this.spanningTree = T;

        this.vertexCount = graph.nv;
        this.edgeCount = graph.ne;

        this.vertexWeights = getVertexWeights();
        this.congestionTree = getCongestionTree();
    }

    // For a graph and its tree, gets the lca values for every edge of the graph
    private int[] lcaResult() {
        TarjanLCA LcaSolver = new TarjanLCA(spanningTree);

        int Index = 0;
        int[] LeftNeighbors = new int[edgeCount];
        int[] RightNeighbors = new int[edgeCount];

        for (int i = 0; i < vertexCount; i++) {
            int NeighborCount = graph.nbrs[i].length;

            for (int j = 0; j < NeighborCount; j++)
                if (i < graph.nbrs[i][j]) { // Check only for ordered edges
                    LeftNeighbors[Index] = i;
                    RightNeighbors[Index] = graph.nbrs[i][j];
                    Index++;
                }
        }

        return LcaSolver.solve(LeftNeighbors, RightNeighbors);
    }

    // Instead of incrementing each vertex of the tree by the value of the path (a,b) in the original graph,
    // will increment vertex a by weight(a,b), increment vertex b by weight(a,b) and decrement vertex lca(a,b)
    // by 2 * weight(a,b). Special cases when lca(a,b) = a or b are treated surprisingly smooth.
    // This will help us have a complexity of O(M) on this segment of the code.

    // Assigns values to vertexWeights
    private double[] getVertexWeights() {
        int[] lcaEdgeValues = lcaResult();
        double[] vertexWeights = new double[vertexCount];

        int Index = 0;

        for (int i = 0; i < vertexCount; i++) {
            for (int j = 0; j < graph.nbrs[i].length; j++)
                if (i < graph.nbrs[i][j]) { // Do only for ordered pairs
                    // A (u,v) edge with LCA and Cost

                    int u = i;
                    int v = graph.nbrs[i][j];
                    int LCA = lcaEdgeValues[Index];
                    double Cost = graph.weights[i][j];

                    vertexWeights[u] += Cost;
                    vertexWeights[v] += Cost;
                    vertexWeights[LCA] -= 2 * Cost;

                    Index++;
                }
        }

        return vertexWeights;
    }

    // Gets the congestionTree
    private Tree getCongestionTree() {
        Tree answer = spanningTree;

        int[] bfsOrdering = TreeUtils.bfsOrder(spanningTree);   // Get the BFS ordering

        for (int i = vertexCount - 1; i > 0; i--) {

            int child = bfsOrdering[i];
            int parent = spanningTree.nodes[child].parent;

            double congestionCost = vertexWeights[child] / spanningTree.nodes[child].length;
            answer.nodes[child].length = congestionCost;

            vertexWeights[parent] += vertexWeights[child];
        }

        return answer;
    }
}
