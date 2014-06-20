/**
 * @file GraphUtils.java
 * @author Alex Reinking <alexander.reinking@gmail.com>
 * @date Thu Jun 5 2014
 *
 * Various static methods to compute useful things on graphs.
 */

package lapsolver.util;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class GraphUtils {
    public static Graph permuteGraph(Graph g, int[] perm) {
        int N = perm.length;
        int M = g.ne;

        int[] posInPermutation = new int[N];
        for (int i = 0; i < N; i++)
            posInPermutation[perm[i]] = i;

        int index = 0;
        int[] src = new int[M];
        int[] dst = new int[M];
        double[] weight = new double[M];

        for (int u = 0; u < N; u++)
            for (int i = 0; i < g.deg[u]; i++) {
                int v = g.nbrs[u][i];
                if (u < v) {
                    src[index] = posInPermutation[u];
                    dst[index] = posInPermutation[v];
                    weight[index] = g.weights[u][i];
                    index++;
                }
            }

        return new Graph(src, dst, weight);
    }

    // given a graph structure that represents a tree, turn it into a rooted tree
    public static Tree toTree(Graph g, int root) {
        // DFS state
        boolean[] visited = new boolean[g.nv];
        int[] stack = new int[g.nv];
        int stack_pos = 1;

        // parent array, to be built
        int[] parent = new int[g.nv];
        double[] weight = new double[g.nv];

        // start at root
        stack[0] = root;
        parent[root] = root;

        // perform DFS
        while(stack_pos > 0) {
            int v = stack[--stack_pos]; // pop vertex
            visited[v] = true;

            // explore children
            for (int i = 0; i < g.deg[v]; i++) {
                int child = g.nbrs[v][i];
                if (visited[child]) continue;
                stack[stack_pos++] = child;
                parent[child] = v;
                weight[child] = g.weights[v][i];
            }
        }

        // build tree from parent array
        return new Tree(parent, weight);
    }

    // build a rooted tree from a graph starting at 0, for convenience
    public static Tree toTree(Graph g) {
        return toTree(g, 0);
    }

    // prints the contents of g to stdout, as an adjacency list
    public static void dump(Graph g) {
        for (int x = 0; x < g.nv; x++) {
            System.out.print(x + " : ");
            for (int i = 0; i < g.deg[x]; i++)
                System.out.print(g.nbrs[x][i] + " ");
            System.out.println();
        }
    }

    /**
     * Compute connected components using a bfs approach.
     * Returns a the characteristic vector of the components,
     * starting to auxiliarySize with 1
     */
    public int[] getComponents(Graph g) {
        int[] order = new int[g.nv];
        int[] comp = new int[g.nv];

        for (int i = 0; i < g.nv; i++) {
            comp[i] = 0;
        }

        int c = 0;
        for (int x = 0; x < g.nv; x++) {
            if (comp[x] == 0) {
                c = c + 1;
                comp[x] = c;
                int ptr = 0;
                int orderLen = 1;
                order[ptr] = x;

                while (ptr < orderLen) {
                    int curNode = order[ptr];
                    for (int i = 0; i < g.deg[curNode]; i++) {
                        int nbr = g.nbrs[curNode][i];
                        if (comp[nbr] == 0) {
                            comp[nbr] = c;
                            order[orderLen++] = nbr;
                        }
                    }
                    ptr++;
                }
            }
        }
        return comp;
    }

    // get a graph from an edge list (for structure) and Graph (for weights)
    // used by low-stretch spanning trees that modify edge weights
    public static Graph getSubgraphStructure(EdgeList edges, Graph graph) {
        Graph graphFromEdges = new Graph(edges);

        // copy weights from old graph
        int[] adjIndex = new int[graph.nv];
        for (int u = 0; u < graph.nv; u++) {
            for (int i = 0; i < graph.deg[u]; i++) {
                adjIndex[graph.nbrs[u][i]] = i;
            }
            for (int i = 0; i < graphFromEdges.deg[u]; i++) {
                int v = graphFromEdges.nbrs[u][i];
                double w = graph.weights[u][adjIndex[v]];
                graphFromEdges.weights[u][i] = w;
            }
        }

        return graphFromEdges;
    }

    // exorcises duplicate edges and self-loops from an edge list
    public static EdgeList sanitizeEdgeList (EdgeList edges) {
        // make scratch space
        final EdgeList auxEdges = new EdgeList(edges);
        ArrayList<Integer> edgeOrder = new ArrayList<>();
        ArrayList<Integer> edgesToKeep = new ArrayList<>();

        // canonize edges so u[i] < v[i]
        for (int i = 0; i < auxEdges.ne; i++) {
            int lo = Math.min(auxEdges.u[i], auxEdges.v[i]);
            int hi = Math.max(auxEdges.u[i], auxEdges.v[i]);
            auxEdges.u[i] = lo;
            auxEdges.v[i] = hi;
        }

        // sort edges so equivalence classes are adjacent
        for (int i = 0; i < auxEdges.ne; i++) {
            edgeOrder.add(i);
        }
        Collections.sort(edgeOrder, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                int cmp = Integer.compare(auxEdges.u[o1], auxEdges.u[o2]);
                if (cmp != 0) return cmp;
                return Integer.compare(auxEdges.v[o1], auxEdges.v[o2]);
            }
        });

        // combine parallel edges
        int firstEqualEdge = -1;
        for (int i = 0; i < auxEdges.ne; i++) {
            int curIndex = edgeOrder.get(i);
            if (auxEdges.u[curIndex] == auxEdges.v[curIndex]) continue;
            if (firstEqualEdge == -1
                    || auxEdges.u[edgeOrder.get(i-1)] != auxEdges.u[curIndex]
                    || auxEdges.v[edgeOrder.get(i-1)] != auxEdges.v[curIndex]) {
                // first time seeing this edge
                edgesToKeep.add(curIndex);
                firstEqualEdge = curIndex;
            }
            else {
                // parallel edge
                auxEdges.weight[firstEqualEdge] = auxEdges.weight[firstEqualEdge] + auxEdges.weight[curIndex];
            }
        }

        // keep selected edges
        EdgeList keptEdges = new EdgeList(edgesToKeep.size());
        int index = 0;
        for (int i : edgesToKeep) {
            keptEdges.u[index] = auxEdges.u[i];
            keptEdges.v[index] = auxEdges.v[i];
            keptEdges.weight[index] = auxEdges.weight[i];
            index++;
        }

        return keptEdges;
    }
}
