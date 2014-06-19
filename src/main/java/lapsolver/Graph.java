/**
 * @file Graph.java
 * @author Daniel Spielman <spielman@cs.yale.edu>
 * @author Cyril Zhang <cyril.zhang@yale.edu>
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
 * graph = Graph(i-1,j-1,v)
 *
 * Elementary methods: dfs, components, bfsWalk, 
 *
 */

package lapsolver;

public class Graph {
    public int nv; // number of vertices
    public int ne; // number of edges

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

    //-----------------------------
    //  backInd[x][i] = the position of x in the neighbor list of u = graph.nbrs[x][i]
    //  nbrs[nbrs[x][i],backInd[x][i]] = x
    public int[][] backInd;

    public Graph() { }

    // constructor from edge list
    public Graph(EdgeList edges) {
        this(edges.u, edges.v, edges.weight);
    }

    // constructor from edge data
    public Graph(int[] src, int[] dst, double[] weight) {
        ne = src.length; // the length of i

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
                throw new Error("Self-loops are not allowed (" + i + "/" + ne + ").");
            }
        }

        // compute max node auxiliarySize
        nv = 0;
        for (int i = 0; i < ne; i++) {
            if (src[i] > nv)
                nv = src[i];
            if (dst[i] > nv)
                nv = dst[i];
        }
        nv++; // largest 0-based auxiliarySize is nv-1

        //  count how many times each node occurs
        deg = new int[nv];
        for (int i = 0; i < nv; i++) {
            deg[i] = 0;
        }

        for (int i = 0; i < ne; i++) {
            deg[src[i]]++;
            deg[dst[i]]++;
        }

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

        for (int i = 0; i < ne; i++) {
            if (src[i] > dst[i]) {
                weights[src[i]][tmpdeg[src[i]]] = weight[i];
                weights[dst[i]][tmpdeg[dst[i]]] = weight[i];

                backInd[src[i]][tmpdeg[src[i]]] = tmpdeg[dst[i]];
                backInd[dst[i]][tmpdeg[dst[i]]] = tmpdeg[src[i]];

                nbrs[src[i]][tmpdeg[src[i]]++] = dst[i];
                nbrs[dst[i]][tmpdeg[dst[i]]++] = src[i];
            }
        }
    }

    // get graph from tree
    public Graph(Tree tree) {
        this(new EdgeList(tree));
    }

    // copy constructor (perform a deep copy)
    public Graph(Graph other) {
        nv = other.nv;
        ne = other.ne;
        deg = other.deg.clone();

        nbrs = new int[nv][];
        weights = new double[nv][];
        backInd = new int[nv][];

        for (int i = 0; i < nv; i++) {
            nbrs[i] = other.nbrs[i].clone();
            weights[i] = other.weights[i].clone();
            backInd[i] = other.backInd[i].clone();
        }
    }
}
