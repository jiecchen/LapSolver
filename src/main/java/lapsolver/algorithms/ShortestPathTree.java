/**
 * @file ShortestPathTree.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Tue Jun 10 2014
 *
 * Slow implementation of Dijkstra's algorithm.
 * Correct, but the boxing/unboxing is a major performance hit.
 *
 */
package lapsolver.algorithms;

import lapsolver.Graph;
import lapsolver.Tree;

import java.util.Comparator;
import java.util.PriorityQueue;


public class ShortestPathTree {
    private final double [] dist;
    private final int [] parent;

    public double[] getDist() {
        return dist;
    }

    public int[] getParent() {
        return parent;
    }

    /**
     * Construct a shortest-path tree for a given graph and source node
     * @param G the input graph
     * @param source the starting node
     */
    public ShortestPathTree(Graph G, int source) {
        dist = new double[G.nv];
        parent = new int[G.nv];

        PriorityQueue<Integer> nextNodes = new PriorityQueue<>(G.nv, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return Double.compare(dist[o1], dist[o2]);
            }
        });

        dist[source] = 0;
        for (int i = 0; i < G.nv; i++) {
            if(i != source) {
                dist[i] = Double.POSITIVE_INFINITY;
                parent[i] = -1;
            }
            nextNodes.add(i);
        }

        boolean[] settled = new boolean[G.nv];
        while (!nextNodes.isEmpty()) {
            int u = nextNodes.poll();
            for (int i = 0; i < G.deg[u]; i++) {
                int v = G.nbrs[u][i];
                if(!settled[v]) {
                    double alt = dist[u] + G.weights[u][i];
                    if(alt < dist[v]) {
                        nextNodes.remove(v);
                        dist[v] = alt;
                        parent[v] = u;
                        nextNodes.add(v);
                    }
                }
            }
            settled[u] = true;
        }

        parent[source] = source; // Follow parent array convention
    }

    /**
     *
     * @return the shortest path tree computed by the constructor
     */
    public Tree getTree() {
        return new Tree(parent, dist);
    }
}
