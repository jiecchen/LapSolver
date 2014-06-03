/**
 * @file UnionFind.java
 * @author Daniel Spielman <spielman@cs.yale.edu>
 * @date Fri May 30 2014
 *
 * A crude Union Find data structure, implemented using linked lists for each set.
 * Should eventually be replaced by a pointer-based structure with path compression.
 *
 */

package lapsolver.algorithms;

public class UnionFind {
    //-----------------------
    //   the number of verts
    //
    public int nv;

    //-----------------------
    //
    // comp[i] is the name of the component of vertex i
    public int[] comp;

    //-----------------------
    //
    // ptr[i] is a pointer to the next element in the same set as item i
    // -1 if it is the end of a list
    public int[] ptr;

    //-----------------------
    //
    // lastPtr[i] is the last item in comp[i], 
    // so ptr[lastPtr[i]] = -1;
    public int[] lastPtr;


    //-----------------------
    //
    public int[] compSize;

    // Call with number of items in the sets
    public UnionFind(int n) {
        this.nv = n;

        comp = new int[n];
        ptr = new int[n];
        lastPtr = new int[n];
        compSize = new int[n];

        for (int i = 0; i < n; i++) {
            ptr[i] = -1;
            compSize[i] = 1;
            comp[i] = i;
            lastPtr[i] = i;
        }
    }

    public int find(int v) {
        return comp[v];
    }

    // merge the components with items u and v
    public void union(int u, int v) {
        int uc = comp[u];
        int vc = comp[v];

        // swap to make vc the smaller one
        if (compSize[vc] > compSize[uc]) {
            int tmp = vc;
            vc = uc;
            uc = tmp;
        }

        // go through the list for vc, re-naming comp for all
        int i = vc;
        while (i != -1) {
            comp[i] = uc;
            i = ptr[i];
        }

        compSize[uc] += compSize[vc];
        compSize[vc] = 0;

        ptr[lastPtr[uc]] = vc;
        lastPtr[uc] = lastPtr[vc];
    }
}
