/**
 * @file ShortestPathTree.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Tue Jun 10 2014
 *
 * Implementation of Dijkstra's algorithm.
 * Correct, but the boxing/unboxing is a major performance hit.
 *
 */
package lapsolver.algorithms;

import lapsolver.Graph;
import lapsolver.Tree;

import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.TreeSet;


public class ShortestPathTree {
    private final double[] dist;
    private final int[] parent;
    private final int[] parentIndex; // index i such that G.nbrs[u][i] = parent[u]
    private final double[] parentWeight;
    private final double radius;

    public double[] getDist() {
        return dist;
    }

    public int[] getParent() {
        return parent;
    }

    public int[] getParentIndex() {
        return parentIndex;
    }

    public double getRadius() { return radius; }

    /**
     * @param G The input graph (weights are resistances).
     * @param source The source (root) of the shortest path tree.
     * @param ignore A vector to ignore vertices (set to non-(-1) to ignore)
     */
    public ShortestPathTree(Graph G, int source, int[] ignore) {
        dist = new double[G.nv];
        parent = new int[G.nv];
        parentIndex = new int[G.nv];
        parentWeight = new double[G.nv];
        double radius = 0.0;

        TreeSet<Integer> nextNodes = new TreeSet<>(new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
            if (dist[o1] == dist[o2]) return Integer.compare(o1, o2);
            return Double.compare(dist[o1], dist[o2]);
            }
        });

        Arrays.fill(dist, Double.POSITIVE_INFINITY);
        Arrays.fill(parent, -1);
        dist[source] = 0;
        parent[source] = source;

        for (int i = 0; i < G.nv; i++) {
            nextNodes.add(i);
        }

        boolean[] settled = new boolean[G.nv];
        while (!nextNodes.isEmpty()) {
            int u = nextNodes.pollFirst();
            if (ignore == null || ignore[u] == -1) {
                for (int i = 0; i < G.deg[u]; i++) {
                    int v = G.nbrs[u][i];
                    if (!settled[v] && (ignore == null || ignore[v] == -1)) {
                        double alt = dist[u] + G.weights[u][i];
                        if (alt < dist[v]) {
                            nextNodes.remove(v);
                            dist[v] = alt;
                            parentWeight[v] = G.weights[u][i];
                            parentIndex[v] = G.backInd[u][i];
                            parent[v] = u;
                            nextNodes.add(v);
                        }
                    }
                }
            }

            if (dist[u] > radius) radius = dist[u];
            settled[u] = true;
        }

        this.radius = radius;
    }

    /**
     * Construct a shortest-path tree for a given graph and source node
     *
     * @param G      the input graph (weights are resistances).
     * @param source the starting node
     */
    public ShortestPathTree(Graph G, int source) {
        this(G, source, null);
    }

    /**
     * @build the shortest path tree using the parent-weight array
     * computed by the constructor
     */
    public Tree getTree() {
        return new Tree(parent, parentWeight);
    }
}
