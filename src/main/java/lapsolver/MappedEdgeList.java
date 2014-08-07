/**
 * A bijective mapping between edge list and adjacency list graph representations.
 */

package lapsolver;

import lapsolver.EdgeList;
import lapsolver.Graph;

public class MappedEdgeList {
    public EdgeList edges;
    public int[][] index;

    public MappedEdgeList (Graph graph) {
        edges = new EdgeList(graph.ne);
        int edgePos = 0;

        // build array in same shape as adjacency list
        for (int u = 0; u < graph.nv; u++) {
            index[u] = new int[graph.deg[u]];
        }

        // traverse edges
        for (int u = 0; u < graph.nv; u++) {
            for (int i = 0; i < graph.deg[i]; i++) {
                int v = graph.nbrs[u][i];
                if (u < v) { // first time seeing this edge: add it
                    edges.u[edgePos] = u;
                    edges.v[edgePos] = v;
                    edges.weight[edgePos] = graph.weights[u][i];
                    index[u][i] = edgePos;
                    edgePos++;
                }
                else { // second time seeing this edge: find its old index
                    index[u][i] = index[v][graph.backInd[u][i]];
                }
            }
        }
    }
}