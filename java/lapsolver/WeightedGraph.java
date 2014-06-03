/**
 * @file WeightedGraph.java
 * @author Daniel Spielman <spielman@cs.yale.edu>
 * @date Sun Dec 22 2013
 *
 * (c) Daniel Spielman, 2013, part of the YINSmat package
 *
 * A sparse weighted static graph class.
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
 * graph = WeightedGraph(i-1,j-1,v)
 *
 * Elementary methods: dfs, components, bfsWalk, 
 *
 */

package lapsolver;

public class WeightedGraph {
    //-----------------------
    //   the number of verts
    //
    public int nv;

    //-----------------------
    //   the number of edges
    //
    public int ne;

    //---------------------
    //  the main graph data,
    //  nbrs[x][i] is ith nbr 
    //       of node x
    //  weights[x][i] is the weight of this edge
    //     (note each appears twice)
    //
    //  deg[x] is the degree of node x

    public int[][] nbrs;
    public double[][] weights;
    public int[] deg;


    // only computed if call compWtedDeg
    public double[] weightedDeg;

    // the sum of all weighted degrees
    public double volume;

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

    public WeightedGraph() { }
    public WeightedGraph(int[] src, int[] dst, double[] weight) {
        fromEdgeList(src, dst, weight);
    }

    /**
     * expects input as
     * [i,j,v] = find(tril(T)),
     * so, each entry of i should be larger than corresp entry of j
     */
    public void fromEdgeList(int[] src, int[] dst, double[] weight) {
        ne = src.length; // the length of i

        dfs = null;
        comp = null;
        seen = null;

        if ((dst.length != ne) || (weight.length != ne)) {
            throw new Error("inputs must all have the same length");
        }

        // Ensure src[i] > dst[i], src[i] != dst[i] \forall i
        for (int i = 0; i < ne; i++) {
            int lower = Math.min(src[i], dst[i]);
            int upper = Math.max(src[i], dst[i]);

            dst[i] = lower;
            src[i] = upper;

            if (src[i] == dst[i]) {
                throw new Error("Self-loops are not allowed.");
            }
        }

        //-------------------------
        // compute max node index
        //
        nv = 0;
        for (int i = 0; i < ne; i++) {
            if (src[i] > nv)
                nv = src[i];
            if (dst[i] > nv)
                nv = dst[i];
        }
        nv++; // largest 0-based index is nv-1

        //-----------------------------------------
        //  count how many times each node occurs
        deg = new int[nv];
        for (int i = 0; i < nv; i++) {
            deg[i] = 0;
        }

        for (int i = 0; i < ne; i++) {
            deg[src[i]]++;
            deg[dst[i]]++;
        }

        //---------------------------
        // build the graph
        int[] tmpdeg = new int[nv];

        nbrs = new int[nv][];
        backInd = new int[nv][];
        weights = new double[nv][];

        for (int i = 0; i < nv; i++) {
            nbrs[i] = new int[deg[i]];
            backInd[i] = new int[deg[i]];
            weights[i] = new double[deg[i]];
            tmpdeg[i] = 0;
        }

        volume = 0;

        for (int i = 0; i < ne; i++) {
            if (src[i] > dst[i]) {
                weights[src[i]][tmpdeg[src[i]]] = weight[i];
                weights[dst[i]][tmpdeg[dst[i]]] = weight[i];
                volume += weight[i];

                backInd[src[i]][tmpdeg[src[i]]] = tmpdeg[dst[i]];
                backInd[dst[i]][tmpdeg[dst[i]]] = tmpdeg[src[i]];

                nbrs[src[i]][tmpdeg[src[i]]++] = dst[i];
                nbrs[dst[i]][tmpdeg[dst[i]]++] = src[i];
            }
        }

        volume *= 2; // double count
    }

    /**
     * set up the graph so that it can be input through fromEdgeList
     * is meant for calls from java: zero indexed, and no clear ordering on i and j
     */
    public void fromMatlab(int[] src, int[] dst, double[] weight) {
        //-------------------------
        // downshift i and j by 1
        //
        for(int i = 0; i < src.length; i++) {
            src[i]--;
            dst[i]--;
        }

        fromEdgeList(src, dst, weight);
    }

    /*
     * The back indices satisfy
     * nbrs[nbrs[x][i],backInd[x][i]] = x
     *
     */
    public void makeBackEdges() {
        backInd = new int[nv][];

        int[] count = new int[nv];

        for (int a = 0; a < nv; a++) {
            backInd[a] = new int[deg[a]];
            count[a] = 0;
        }

        for (int a = 0; a < nv; a++) {
            for (int i = 0; i < deg[a]; i++) {
                int nbr = nbrs[a][i];
                backInd[nbr][count[nbr]++] = i;
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
     * starting to index with 1
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
       this turns it into a Tree class object via a parent array
    */
    public Tree treeToTree() {
        int[] pArray = treeToArray();

        // now, compute all the edge weights
        double[] wt = new double[nv];

        for (int u = 0; u < nv; u++) {
            int v = pArray[u];
            if (v != u) {
                for (int j = 0; j < deg[v]; j++)
                    if (nbrs[v][j] == u)
                        wt[u] = weights[v][j];
            }
        }
        return new Tree(pArray, wt);
    }

    public void dump() {
        for (int x = 0; x < nv; x++) {
            System.out.print(x + " : ");
            for (int i = 0; i < deg[x]; i++)
                System.out.print(nbrs[x][i] + " ");
            System.out.println();
        }

    }

    public void compWtedDeg() {
        weightedDeg = new double[nv];
        for (int x = 0; x < nv; x++) {
            double sum = 0;
            for (int i = 0; i < deg[x]; i++)
                sum += weights[x][i];
            weightedDeg[x] = sum;
        }
    }

    /*
     * Return a 3-by-m ijv vector.
     * In Matlab, we can turn this into an adjacency matrix by
     * matrix = sparse(i,j,v,n,n); matrix = matrix + matrix';
     *
     * so, this just spits out the lower-triangular part.
     * For the full matrix, use toIJVsym
     */
    public double[][] toIJV() {
        double[] i = new double[ne];
        double[] j = new double[ne];
        double[] v = new double[ne];

        int ptr = 0;

        for (int x = 0; x < nv; x++) {
            for (int y = 0; y < deg[x]; y++) {
                if (nbrs[x][y] < x) {
                    i[ptr] = (double) x;
                    j[ptr] = (double) nbrs[x][y];
                    v[ptr] = weights[x][y];
                    ptr++;
                }
            }
        }

        double[][] ijv = new double[3][];
        ijv[0] = i;
        ijv[1] = j;
        ijv[2] = v;

        return ijv;
    }

    /*
     * Return a 3-by-m ijv vector.
     * In Matlab, we can turn this into an adjacency matrix by
     * matrix = sparse(ijv(1,:)+1,ijv(2,:)+1,v,n,n);
     *
     * so, this just spits out the lower-triangular part.
     * For the full matrix, use toIJVsym
     */
    public double[][] toIJVsym() {
        double[] i = new double[2 * ne];
        double[] j = new double[2 * ne];
        double[] v = new double[2 * ne];

        int ptr = 0;

        for (int x = 0; x < nv; x++) {
            for (int y = 0; y < deg[x]; y++) {
                i[ptr] = (double) x;
                j[ptr] = (double) nbrs[x][y];
                v[ptr] = weights[x][y];
                ptr++;
            }
        }

        double[][] ijv = new double[3][];
        ijv[0] = i;
        ijv[1] = j;
        ijv[2] = v;

        return ijv;
    }

    /*
     * Produces a copy of this graph
     */
    public WeightedGraph copy() {
        WeightedGraph G = new WeightedGraph();

        G.nv = nv;
        G.ne = ne;
        G.deg = new int[nv];

        System.arraycopy(deg, 0, G.deg, 0, nv);

        G.nbrs = new int[nv][];
        G.weights = new double[nv][];

        for (int a = 0; a < nv; a++) {
            G.nbrs[a] = new int[deg[a]];
            G.weights[a] = new double[deg[a]];
            for (int i = 0; i < deg[a]; i++) {
                G.nbrs[a][i] = nbrs[a][i];
                G.weights[a][i] = weights[a][i];
            }
        }

        return G;
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
