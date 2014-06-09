/**
 * @file GraphUtils.java
 * @author Alex Reinking <alexander.reinking@gmail.com>
 * @date Thu Jun 5 2014
 *
 * Various static methods to compute useful things on graphs.
 */

package lapsolver.util;

import lapsolver.Graph;
import lapsolver.Tree;

public class GraphUtils {
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
     * starting to index with 1
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
}
