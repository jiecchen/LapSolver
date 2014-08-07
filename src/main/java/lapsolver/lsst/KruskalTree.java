/**
 * @file KruskalTree.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Tue Jun 3 2014
 *
 * A not-so-low-stretch spanning tree strategy, which just uses Kruskal's algorithm.
 * This will not produce good low stretch spanning trees.
 *
 */
package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.UnionFind;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class KruskalTree implements SpanningTreeStrategy {
    @Override
    public Tree getTree(Graph graph) {
        return new Tree(getTreeEdges(graph));
    }

    public EdgeList getTreeEdges(Graph graph) {
        EdgeList edges = new EdgeList(graph);
        return UnionFindTreeBuilder.buildTree(edges, graph.nv, edges.weight);
    }
}
