package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.Stretch;

import java.util.LinkedList;
import java.util.Queue;

/**
 * Created by serbanstan on 8/4/14.
 * Starts bfs searches in each vertex and returns the tree with minimum stretch. Should only work for graphs with
 * weight 1 edges.
 */

public class BFSTree implements SpanningTreeStrategy {
    @Override
    public Tree getTree(Graph graph) {
        double bestStretch = -1;
        EdgeList finalTree = new EdgeList(graph.nv - 1);

        for (int start = 0; start < graph.nv; start++) {
            int currentEdge = 0;
            EdgeList treeEdges = new EdgeList(graph.nv - 1);

            int[] order = new int[graph.nv];
            int orderPos = 0;
            boolean[] visited = new boolean[graph.nv];

            Queue<Integer> bfsQueue = new LinkedList<>();
            bfsQueue.add(start);

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

            Tree tree = new Tree(treeEdges);
            Stretch.StretchResult treeStretch = Stretch.compute(graph, tree);

            if (treeStretch.total < bestStretch || bestStretch == -1) {
                bestStretch = treeStretch.total;
                finalTree = new EdgeList(treeEdges);
            }
        }

        return new Tree(finalTree);
    }

}
