/**
 * @file BFSTree.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Aug 6 2014
 *
 * Starts bfs searches in each vertex and returns the tree with minimum stretch. Should only work for graphs with
 * weight 1 edges.
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.Stretch;

import java.util.LinkedList;
import java.util.Queue;

public class BFSTree implements SpanningTreeStrategy {
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
        int currentEdge = 0;
        EdgeList treeEdges = new EdgeList(graph.nv - 1);

        int[] order = new int[graph.nv];
        int orderPos = 0;
        boolean[] visited = new boolean[graph.nv];

        Queue<Integer> bfsQueue = new LinkedList<>();
        bfsQueue.add(root);

        while (!bfsQueue.isEmpty()) {
            int u = bfsQueue.poll();
            visited[u] = true;

            order[orderPos++] = u;

            for (int i = 0; i < graph.deg[u]; i++)
                if (visited[graph.nbrs[u][i]] == false) {
                    bfsQueue.add(graph.nbrs[u][i]);
                    visited[graph.nbrs[u][i]] = true;

                    treeEdges.u[currentEdge] = u;
                    treeEdges.v[currentEdge] = graph.nbrs[u][i];
                    treeEdges.weight[currentEdge] = graph.weights[u][i];
                    currentEdge++;
                }
        }

        return new Tree(treeEdges);
    }

}
