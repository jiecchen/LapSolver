/**
 * @file DijkstraTree.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Aug 6 2014
 *
 * Finds the lowest-stretch shortest path tree.
 */

package lapsolver.lsst;

import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.ShortestPathTree;
import lapsolver.algorithms.Stretch;

public class DijkstraTree implements SpanningTreeStrategy {
    @Override
    public Tree getTree(Graph graph) {
        double bestStretch = 0;
        Tree bestTree = null;

        for (int root = 0; root < graph.nv; root++) {
            Tree candidateTree = getCandidate(graph, root);
            Stretch.StretchResult treeStretch = Stretch.compute(graph, candidateTree);

            if (bestTree == null || treeStretch.total < bestStretch) {
                bestStretch = treeStretch.total;
                bestTree = candidateTree;
            }
        }
        return bestTree;
    }

    public Tree getCandidate(Graph graph, int root) {
        ShortestPathTree spt = new ShortestPathTree(graph, root);
        return spt.getTree();
    }

}
