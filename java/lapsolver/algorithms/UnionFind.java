/**
 * @file UnionFind.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Tue Jun 3 2014
 *
 * A disjoint forest data structure implementing path compression.
 * Union-by-rank is not implemented. The overhead seems not to be worth the asymptotic gain, after testing.
 */

package lapsolver.algorithms;

public class UnionFind {
    //-----------------------
    // the number of vertices
    //
    public int nv;

    //-----------------------
    // a pointer to its parent vertex, or -1 if none
    //
    public int[] ptr;

    // Call with number of items in the sets
    public UnionFind(int n) {
        this.nv = n;

        ptr = new int[n];

        for (int i = 0; i < n; i++) {
            ptr[i] = i;
        }
    }

    // name of component that contains vertex v
    public int find(int v) {
        if (ptr[v] == v) {
            // at root of this tree
            return v;
        } else {
            // follow pointer; do path compression
            ptr[v] = find(ptr[v]);
            return ptr[v];
        }
    }

    // merge the components with vertices u and v
    public void union(int u, int v) {
        int ru = find(u), rv = find(v); // go to roots of u and v
        if(ru != rv) ptr[ru] = rv; // attach if different
    }
}
