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
import java.util.ArrayList;

public class TarjanLCA {
    private Tree tree;

    // state variables for algorithm
    private UnionFind unionFind;
    private boolean[] black;
    private int[] ancestor;

    // store queries as a subgraph
    private int nq;
    private int[][] queries;

    // initialize state
    public TarjanLCA (Tree tree) {
        this.tree = tree;
        unionFind = new UnionFind(tree.nv);
        black = new boolean[tree.nv];
        ancestor = new int[tree.nv];
        queries = new int[tree.nv][];
    }

    // do preprocessing on queries
    public void initQueries (int[] a, int[] b) {
        if (a.length != b.length) {
            throw new Error("LCA query arrays should have same length");
        }
        nq = a.length;

        // one pass to compute degrees
        int[] queryDegree = new int[tree.nv];
        for (int i = 0; i < nq; i++) {
            queryDegree[a[i]]++;
            queryDegree[b[i]]++;
        }

        // init query array sizes
        for (int i = 0; i < tree.nv; i++) {
            queries[i] = new int[queryDegree[i]];
        }

        // populate queries
        for (int i = 0; i < nq; i++) {
            queries[ a[i] ][ --queryDegree[a[i]] ] = b[i];
            queries[ b[i] ][ --queryDegree[b[i]] ] = a[i];
        }
    }

    // solve queries
    public int[] solve() {
        dfs(tree.root);
        return new int[nq];
    }

    // main DFS procedure (from Wikipedia article)
    private void dfs(int u) {
        ancestor[u] = u;
        for (int v : tree.nodes[u].kids) {
            dfs(v);
            unionFind.union(u, v);
            ancestor[unionFind.find(u)] = u;
        }
        black[u] = true;

        for (int i = 0; i < queries[u].length; i++) {
            int v = queries[u][i];
            if (black[v]) {
                int lca = ancestor[unionFind.find(v)];
                System.out.println("LCA " + u + " " + v + " = " + lca);
            }
        }
    }

}
