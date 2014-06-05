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
    private int[] ptr;

    // internal state variables for find
    private int[] find_stack;
    private int find_stack_pos;

    // Call with number of items in the sets
    public UnionFind(int n) {
        this.nv = n;

        ptr = new int[n];
        find_stack = new int[n];

        for (int i = 0; i < n; i++) {
            ptr[i] = i;
        }
    }

    // name of component that contains vertex v
    public int find(int v) {
        find_stack_pos = 0;

        // follow pointers until self-loop
        while(ptr[v] != v) {
            find_stack[find_stack_pos++] = v;
            v = ptr[v];
        }

        // path compression
        while(--find_stack_pos > 0) {
            find_stack[find_stack_pos] = v;
        }

        return v;
    }

    // merge the components with vertices u and v
    public void union(int u, int v) {
        int ru = find(u), rv = find(v); // go to roots of u and v
        if(ru != rv) ptr[ru] = rv; // attach if different
    }
}
