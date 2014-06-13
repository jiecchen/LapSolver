/**
 * @file KruskalTree.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @author Serban Stan <serban.stan@yale.edu>
 * @date Tue Jun 3 2014
 *
 * A dummy low-stretch spanning tree strategy, which just uses Kruskal's algorithm. This will not produce
 * good low stretch spanning trees.
 *
 * Still a stub.
 *
 */
package lapsolver.lsst;

import com.google.common.collect.ContiguousSet;
import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.UnionFind;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class KruskalTree implements SpanningTreeStrategy {
    @Override
    public Tree getTree(Graph graph) {
        UnionFind disjointSet = new UnionFind(graph.nv);
        EdgeList treeEdges = new EdgeList(graph.nv-1);

        // Compute index
        EdgeList inputEdges = new EdgeList(graph);
        final double[] weights = inputEdges.weight;
        List<Integer> indices = ContiguousSet.create(
                Range.closedOpen(0, weights.length),
                DiscreteDomain.integers()).asList();

        Collections.sort(indices, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return Double.compare(weights[o1], weights[o2]);
            }
        });

        int currentEdge = 0;
        for (Integer edge : indices) {
            int u = inputEdges.u[edge];
            int v = inputEdges.v[edge];

            if(disjointSet.find(u) != disjointSet.find(v)) {
                disjointSet.union(u,v);

                treeEdges.u[currentEdge] = u;
                treeEdges.v[currentEdge] = v;
                treeEdges.weight[currentEdge++] = weights[edge];
            }
        }
        return new Tree(treeEdges);
    }
}
