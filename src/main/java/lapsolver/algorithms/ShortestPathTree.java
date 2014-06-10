/**
 * @file ShortestPathTree.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Tue Jun 10 2014
 *
 * Crappy implementation of Dijkstra's algorithm. Probably has bugs, but is pretty close.
 *
 */
package lapsolver.algorithms;

import lapsolver.Graph;
import lapsolver.Tree;

import java.util.Comparator;
import java.util.HashSet;
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

    public ShortestPathTree(Graph input, int source) {
        dist = new double[input.nv];
        parent = new int[input.nv];

        PriorityQueue<Integer> nextNodes = new PriorityQueue<>(input.nv, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return (int) Math.signum(dist[o1] - dist[o2]);
            }
        });

        dist[source] = 0;
        for (int i = 0; i < input.nv; i++) {
            if(i != source) {
                dist[i] = Double.POSITIVE_INFINITY;
                parent[i] = -1;
            }
            nextNodes.add(i);
        }

        HashSet<Integer> settled = new HashSet<>();
        while (!nextNodes.isEmpty()) {
            int u = nextNodes.poll();
            for (int v : input.nbrs[u]) {
                if(!settled.contains(v)) {
                    double alt = dist[u] + input.weights[u][v];
                    if(alt < dist[v]) {
                        nextNodes.remove(v);
                        dist[v] = alt;
                        parent[v] = u;
                        nextNodes.add(v);
                    }
                }
            }
            settled.add(u);
        }
        parent[source] = source; // Follow parent array convention
    }

    public Tree getTree() {
        return new Tree(parent);
    }
}
