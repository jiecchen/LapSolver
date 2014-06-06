/**
 * @file Congestion.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @date Wed Jun 4 2014
 *
 * For an input graph G and tree T, returns the tree with congestion values on each edge (.getCongestionTree())
 * Complexity: O( M log N )
 *
 * Instead of incrementing each vertex of the tree by the value of the path (a,b) in the original graph,
 * will increment vertex a by weight(a,b), increment vertex b by weight(a,b) and decrement vertex lca(a,b)
 * by 2 * weight(a,b). Special cases when lca(a,b) = a or b are treated surprisingly smooth.
 * This will help us have a complexity of O(M) on this segment of the code.
 */

package lapsolver.algorithms;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.util.TreeUtils;

public class Congestion {
    public static Tree compute(Graph graph, Tree spanningTree) {
        double[] vertexWeights = computeVertexWeights(spanningTree, graph);
        return computeCongestionTree(spanningTree, vertexWeights);
    }

    // Assigns values to vertexWeights
    private static double[] computeVertexWeights(Tree spanningTree, Graph graph) {
        // Compute the LCAs
        TarjanLCA lcaSolver = new TarjanLCA(spanningTree);
        EdgeList edges = new EdgeList(graph);
        int[] lcaEdgeValues = lcaSolver.solve(edges.u, edges.v);

        double[] vertexWeights = new double[graph.nv];

        int index = 0;

        for (int i = 0; i < graph.nv; i++) {
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
    private static Tree computeCongestionTree(Tree spanningTree, double[] vertexWeights) {
        Tree answer = new Tree(spanningTree);

        int[] bfsOrdering = TreeUtils.bfsOrder(spanningTree);   // Get the BFS ordering

        for (int i = spanningTree.nv - 1; i > 0; i--) {
            int child = bfsOrdering[i];
            int parent = spanningTree.getNode(child).getParent().getId();

            // congestion cost
            answer.getNode(child).setLength(vertexWeights[child] / spanningTree.getNode(child).getLength());

            vertexWeights[parent] += vertexWeights[child];
        }

        return answer;
    }
}
