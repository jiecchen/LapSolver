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
import lapsolver.EdgeList;
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

        this.vertexWeights = getVertexWeights();
        this.congestionTree = getCongestionTree();
    }

    // For a graph and its tree, gets the lca values for every edge of the graph
    private int[] lcaResult() {
        TarjanLCA lcaSolver = new TarjanLCA(spanningTree);
        EdgeList edges = new EdgeList(graph);

        return lcaSolver.solve(edges.u, edges.v);
    }

    // Instead of incrementing each vertex of the tree by the value of the path (a,b) in the original graph,
    // will increment vertex a by weight(a,b), increment vertex b by weight(a,b) and decrement vertex lca(a,b)
    // by 2 * weight(a,b). Special cases when lca(a,b) = a or b are treated surprisingly smooth.
    // This will help us have a complexity of O(M) on this segment of the code.

    // Assigns values to vertexWeights
    private double[] getVertexWeights() {
        int[] lcaEdgeValues = lcaResult();
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

    // Gets the congestion tree
    private Tree getCongestionTree() {
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
