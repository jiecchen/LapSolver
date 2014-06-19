/**
 * @file UnweightedGraph.java
 * @author Daniel Spielman <spielman@cs.yale.edu>
 * @date Sun Dec 22 2013
 *
 * (c) Daniel Spielman, 2013, part of the YINSmat package
 *
 * A sparse UNWEIGHTED static graph class.
 * keeps an adj list for each vertex, 
 * not necessarily sorted
 * produces edge back indices only when needed
 *
 * Meant for use inside matlab, but this causes 0/1 indexing problems!
 * Matlab indexes from 1, this code from 0.
 * Matlab must do all the conversions.
 *
 * To create in matlab, use
 * [i,j,v] = find(tril(T))
 * graph = Graph(i-1,j-1,v)
 *
 * Elementary methods: dfs, components, bfsWalk, 
 *
 */

package lapsolver;

public class UnweightedGraph {

    //-----------------------
    //   the number of verts
    //
    public int nv;

    //---------------------
    //  the main graph data,
    //  nbrs[x][i] is ith nbr 
    //       of node x
    //
    //  deg[x] is the degree of node x

    public int[][] nbrs;
    public int[] deg;


    //-----------------------------
    //  nbrs[nbrs[x][i],backInd[x][i]] = x
    public int[][] backInd;
    //--------------------
    //  a bfs ordering,
    //  and bfs depth
    //
    public int[] bfs;
    public int[] depth;
    public boolean[] seen;
    public int[] comp;
    //--------------------
    //  a dfs ordering
    //  not initialized
    //
    private int[] dfs;
    private int dfsPtr; // a utility
    //------------------------
    //  the parent in the dfs tree
    //  initialized by dfs()
    //
    private int[] parent;

    public UnweightedGraph(int[] i, int[] j) {
        this.setGraph(i, j);
    }

    /**
     * expects input as
     * [src,dst] = find(tril(T)),
     * so, each entry of src should be larger than corresponding entry of dst
     */
    public void setGraph(int [] src, int [] dst) {
        int len = src.length; // the length of i

        dfs = null;
        comp = null;
        seen = null;

        if (dst.length != len) {
            throw new Error("inputs must all have the same length");
        }

        //-------------------------
        // compute max node auxiliarySize
        //
        nv = 0;
        for (int i = 0; i < len; i++) {
            if (src[i] > nv)
                nv = src[i];
            if (dst[i] > nv)
                nv = dst[i];
        }


        //-----------------------------------------------------------
        // downshift i and j by 1 to convert from MATLAB conventions
        //
        for (int i = 0; i < len; i++) {
            src[i] = src[i] - 1;
            dst[i] = dst[i] - 1;

            // report error if self-loops
            if (src[i] == dst[i]) {
                throw new Error("Self-loops are not allowed.");
            }
        }


        //-----------------------------------------
        //  count how many times each node occurs
        //
        deg = new int[nv];
        for (int i = 0; i < nv; i++)
            deg[i] = 0;

        for (int i = 0; i < len; i++) {
            if (src[i] > dst[i]) {
                deg[src[i]]++;
                deg[dst[i]]++;
            }
        }

        //---------------------------
        // build the graph

        int[] tmpdeg = new int[nv];

        /*
         * The back indices satisfy
         * nbrs[nbrs[x][i],backInd[x][i]] = x
         *
         */

        nbrs = new int[nv][];
        backInd = new int[nv][];

        for (int i = 0; i < nv; i++) {
            nbrs[i] = new int[deg[i]];
            backInd[i] = new int[deg[i]];
            tmpdeg[i] = 0;
        }

        for (int i = 0; i < len; i++) {
            if (src[i] > dst[i]) {
                backInd[src[i]][tmpdeg[src[i]]] = tmpdeg[dst[i]];
                backInd[dst[i]][tmpdeg[dst[i]]] = tmpdeg[src[i]];

                nbrs[src[i]][tmpdeg[src[i]]++] = dst[i];
                nbrs[dst[i]][tmpdeg[dst[i]]++] = src[i];
            }
        }
    }

    /**
     * walk in BFS order, setting the field bfs
     * assumes connected
     */
    public int[] bfsWalk(int root) {
        bfs = new int[nv];
        depth = new int[nv];
        seen = new boolean[nv];

        for (int i = 0; i < nv; i++) {
            seen[i] = false;
        }

        depth[root] = 0;
        bfs[0] = root;
        seen[root] = true;

        int bfsPtr = 1;
        int curPtr = 1;
        int curNode = root;

        while (curPtr < nv) {
            int nbr;
            for (int i = 0; i < deg[curNode]; i++) {
                nbr = nbrs[curNode][i];
                if (!seen[nbr]) {
                    bfs[bfsPtr++] = nbr;
                    depth[nbr] = depth[curNode] + 1;
                    seen[nbr] = true;
                }
            }
            curNode = bfs[curPtr++];
        }

        return bfs;
    }


    /**
     * Compute dfs by a recursive algorithm
     * could easily break the stack
     */
    public int[] dfs() {
        int c = 0;

        seen = new boolean[nv];
        comp = new int[nv];
        for (int i = 0; i < nv; i++) {
            seen[i] = false;
            comp[i] = 0;
        }

        dfsPtr = 0;
        dfs = new int[nv];
        parent = new int[nv];

        for (int x = 0; x < nv; x++) {
            if (!seen[x]) {
                dfsSub(x, c);
                c = c + 1;
            }
        }

        return dfs;
    }

    public int[] getDfs() {
        return dfs;
    }

    /**
     * Compute connected components using a bfs approach.
     * Returns a the characteristic vector of the components,
     * starting to auxiliarySize with 1
     */
    public int[] components() {

        int[] order = new int[nv];

        seen = new boolean[nv];
        comp = new int[nv];

        for (int i = 0; i < nv; i++) {
            comp[i] = 0;
        }


        int c = 0;
        for (int x = 0; x < nv; x++) {
            if (comp[x] == 0) {
                c = c + 1;
                comp[x] = c;
                int ptr = 0;
                int orderLen = 1;
                order[ptr] = x;

                while (ptr < orderLen) {
                    int curNode = order[ptr];
                    for (int i = 0; i < deg[curNode]; i++) {
                        int nbr = nbrs[curNode][i];
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

    public int[] treeToArray() {

        dfs();

        int[] pArray = new int[nv];

        // set all others to parent
        System.arraycopy(parent, 0, pArray, 0, nv);

        // set root to itself
        pArray[dfs[0]] = dfs[0];

        return pArray;

    }

    
    
    /* if the graph is a tree,
       this turns it into a parent array
    */

    /*
     * Find a maximal indep set by first choosing
     * vertices of low degree
     *
     * return the char vector of the set
     */
    public int[] bigIndepSet() {
        int[] charv = new int[nv];

        // 1 if in mis, 0 if not, -1 if unassigned
        for (int i = 0; i < nv; i++) {
            charv[i] = -1;
        }


        // form link lists for every degree, by ints
        // -1 if empty
        // list[v] is the next item after v
        int[] list = new int[nv];
        for (int i = 0; i < nv; i++) {
            list[i] = -1;
        }

        // also, keep a head of each list
        // -1 if empty
        int[] head = new int[nv];
        for (int i = 0; i < nv; i++) {
            head[i] = -1;
        }

        for (int v = 0; v < nv; v++) {
            int d = deg[v];
            list[v] = head[d];
            head[d] = v;
        }


        // now, start making the mis
        for (int d = 0; d < nv; d++) {
            int v = head[d];
            while (v != -1) {
                if (charv[v] == -1) {
                    charv[v] = 1;
                    for (int i = 0; i < deg[v]; i++) {
                        charv[nbrs[v][i]] = 0;
                    }
                }
                v = list[v];
            }
        }

        return charv;
    }

    public void dump() {

        for (int x = 0; x < nv; x++) {
            System.out.print(x + " : ");
            for (int i = 0; i < deg[x]; i++)
                System.out.print(nbrs[x][i] + " ");
            System.out.println();
        }

    }

    private void dfsSub(int x, int c) {
        dfs[dfsPtr++] = x;
        seen[x] = true;
        comp[x] = c;
        for (int i = 0; i < deg[x]; i++) {
            if (!seen[nbrs[x][i]]) {
                parent[nbrs[x][i]] = x;
                dfsSub(nbrs[x][i], c);
            }
        }
    }
}
