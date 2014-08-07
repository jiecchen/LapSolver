/**
 * @file UnionFindTreeBuilder.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Aug 7 2014
 *
 * A tree strategy that takes an ordering of edges, and adds them in the style of Kruskal's algorithm.
 * Complexity: O(m * alpha(n)) if ordering given
 *             O(m log n) if priorities given
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.algorithms.UnionFind;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class UnionFindTreeBuilder {
    /**
     * Builds a tree from an edge ordering.
     * @param edges A graph in EdgeList form.
     * @param nv The number of vertices in the graph.
     * @param order A permutation of edge indices.
     * @return The tree, as an EdgeList with nv-1 elements.
     */
    public static EdgeList buildTree(EdgeList edges, int nv, int[] order) {
        UnionFind uf = new UnionFind(nv);
        EdgeList treeEdges = new EdgeList(nv-1);
        int edgeIndex = 0;

        for (int e : order) {
            int u = edges.u[e];
            int v = edges.v[e];

            if (uf.find(u) != uf.find(v)) {
                uf.union(u, v);

                treeEdges.u[edgeIndex] = u;
                treeEdges.v[edgeIndex] = v;
                treeEdges.weight[edgeIndex++] = edges.weight[e];
            }
        }

        return treeEdges;
    }

    /**
     * Builds a tree from an edge adding times.
     * @param edges A graph in EdgeList form.
     * @param nv The number of vertices in the graph.
     * @param time Times to add vertices.
     * @return The tree, as an EdgeList with nv-1 elements.
     */
    public static EdgeList buildTree(EdgeList edges, int nv, final double[] time) {
        ArrayList<Integer> order = new ArrayList<>(edges.ne);
        for (int i = 0; i < edges.ne; i++) {
            order.add(i);
        }

        Collections.sort(order, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return Double.compare(time[o1], time[o2]);
            }
        });

        int[] order_prim = new int[edges.ne];
        for (int i = 0; i < edges.ne; i++) {
            order_prim[i] = order.get(i);
        }

        return buildTree(edges, nv, order_prim);
    }
}
