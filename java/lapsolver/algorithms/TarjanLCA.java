/**
 * @file TarjanLCA.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Tue Jun 3 2014
 *
 * An algorithm for computing offline LCA queries in a tree.
 * Complexity: amortized O( q alpha(n) )
 * Reference: http://en.wikipedia.org/wiki/Tarjan's_off-line_lowest_common_ancestors_algorithm
 */

package lapsolver.algorithms;

import lapsolver.Tree;

public class TarjanLCA {
    private Tree tree;

    // state variables for algorithm
    private UnionFind unionFind;
    private boolean[] black;
    private int[] ancestor;

    // store queries as a subgraph
    private int nq;
    private int[][] queries, queryIndices;
    private int[] answer;

    // initialize state
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
        nq = a.length;
        answer = new int[nq];

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

        dfs(tree.root);

        return answer;
    }

    // main DFS procedure (from Wikipedia article)
    private void dfs(int u) {
        ancestor[u] = u;
        for (int v : tree.nodes[u].children) {
            dfs(v);
            unionFind.union(u, v);
            ancestor[unionFind.find(u)] = u;
        }
        black[u] = true;

        for (int i = 0; i < queries[u].length; i++) {
            int v = queries[u][i];
            if (black[v]) {
                answer[ queryIndices[u][i] ] = ancestor[unionFind.find(v)];
            }
        }
    }

}
