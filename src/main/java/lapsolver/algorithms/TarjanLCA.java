/**
 * @file TarjanLCA.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Tue Jun 3 2014
 *
 * An algorithm for computing offline LCA queries in a tree.
 * Complexity: amortized O( (numRemoved + q) alpha(numRemoved) )
 * Reference: http://en.wikipedia.org/wiki/Tarjan's_off-line_lowest_common_ancestors_algorithm
 */

package lapsolver.algorithms;

import lapsolver.Tree;
import lapsolver.util.TreeUtils;

public class TarjanLCA {
    private final Tree tree;

    // state variables for algorithm
    private final UnionFind unionFind;
    private final boolean[] black;
    private final int[] ancestor;

    private final int[][] queries;
    private final int[][] queryIndices;

    /**
     * Initialize the scratch space needed for the algorithm.
     * @param tree The input tree (weights are conductances or resistances).
     */
    public TarjanLCA (Tree tree) {
        this.tree = tree;

        // set up algorithm state
        unionFind = new UnionFind(tree.nv);
        black = new boolean[tree.nv];
        ancestor = new int[tree.nv];

        // set up adjacency lists for queries
        queries = new int[tree.nv][];
        queryIndices = new int[tree.nv][];
    }

    // preprocess, solve, return answers
    public int[] solve (int[] a, int[] b) {
        if (a.length != b.length) {
            throw new Error("LCA query arrays should have same length");
        }
        int nq = a.length;
        int[] answer = new int[nq];

        // one pass to compute degrees
        int[] queryDegree = new int[tree.nv];
        for (int i = 0; i < nq; i++) {
            queryDegree[a[i]]++;
            queryDegree[b[i]]++;
        }

        // init query array sizes
        for (int i = 0; i < tree.nv; i++) {
            queries[i] = new int[queryDegree[i]];
            queryIndices[i] = new int[queryDegree[i]];
        }

        // populate queries
        for (int i = 0; i < nq; i++) {
            queries[ a[i] ][ --queryDegree[a[i]] ] = b[i];
            queries[ b[i] ][ --queryDegree[b[i]] ] = a[i];

            queryIndices[ a[i] ][ queryDegree[a[i]] ] = i;
            queryIndices[ b[i] ][ queryDegree[b[i]] ] = i;
        }

        // do DFS
        for (int i = 0; i < tree.nv; i++) {
            ancestor[i] = i;
        }

        int[] order = TreeUtils.dfsOrder(tree);
        int[] childrenVisited = new int[tree.nv];

        for (int i = tree.nv-1; i >= 0; i--) {
            int v = order[i];
            int parent = tree.parent[v];

            if (v != tree.root) {
                childrenVisited[parent]++;
            }

            if (childrenVisited[v] == tree.children[v].length) {
                black[v] = true;

                for (int j = 0; j < queries[v].length; j++) {
                    int child = queries[v][j];
                    if (black[child]) {
                        answer[ queryIndices[v][j] ] = ancestor[unionFind.find(child)];
                    }
                }

                unionFind.union(parent, v);
                ancestor[unionFind.find(parent)] = parent;
            }
        }

        return answer;
    }

}
