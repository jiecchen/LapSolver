/**
 * @file   WeightedGraph.java
 * @author Daniel Spielman <spielman@cs.yale.edu>
 * @date   Sun Dec 22 2013
 *
 * (c) Daniel Spielman, 2013, part of the YINSmat package
 *
 * A sparse weighted static graph class.
 * keeps an adj list for each vertex, 
 * not nec sorted
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

package yins;

import java.lang.*;
import java.util.*;
import java.io.*;


public class WeightedGraph 
{

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
    public double[] wtedDeg;

    public double volume; // the sum of the weights

    
    //-----------------------------
    //  nbrs[nbrs[x][i],backInd[x][i]] = x
    public int[][] backInd;
    
    //--------------------
    //  a dfs ordering
    //  not initialized
    //
    private int[] dfs;

    private int dfsPtr; // a utility


    //--------------------
    //  a bfs ordering,
    //  and bfs depth
    //
    public int[] bfs;
    public int[] depth;


    //------------------------
    //  the parent in the dfs tree
    //  initialized by dfs()
    //
    private int[] parent;
    private int root;

    public boolean[] seen;
    public int[] comp;


    public WeightedGraph (int[] i, int[] j, double[] v) {
	this.setGraph(i,j,v);
    }

    

    public WeightedGraph () {
	;
    }


    
    /**
     *  expects input as
     *  [i,j,v] = find(tril(T)),
     *  so, each entry of i should be larger than corresp entry of j
     */
    public void setGraph (int[] i, int[] j, double[] v) {
	int len; // the length of i

	dfs = null;
	comp = null;
	seen = null;

	
	len = i.length;
	if ((j.length != len) || (v.length != len)) {
	    throw new Error("inputs must all have the same length");
	}

	//-------------------------
	// compute max node index
	//
	this.nv = 0;
	this.ne = 0;
	for (int a = 0; a < len; a++) {
	    if (i[a] > nv)
		nv = i[a];
	    if (j[a] > nv)
		nv = i[a];
	}


	//-------------------------
	// downshift i and j by 1
	// 
	for (int a = 0; a < len; a++) {
	    i[a] = i[a] - 1;
	    j[a] = j[a] - 1;

	    // report error if self-loops
	    if (i[a] == j[a]) {
		throw new Error("Self-loops are not allowed.");
	    }
	}


	//-----------------------------------------
	//  count how many times each node occurrs
	
	deg = new int[nv];
	for (int a = 0; a < nv; a++) {
	    deg[a] = 0;
	}

	for (int a = 0; a < len; a++) {
	    if (i[a] > j[a]) {
		ne++;
		deg[i[a]] = deg[i[a]] + 1;
		deg[j[a]] = deg[j[a]] + 1;
	    }
	}


	//---------------------------
	// build the graph

	int[] tmpdeg = new int[nv];

	nbrs = new int[nv][];
	weights = new double[nv][];

	for (int a = 0; a < nv; a++) {
	    nbrs[a] = new int[deg[a]];
	    weights[a] = new double[deg[a]];
	    tmpdeg[a] = 0;
	}

	volume = 0;
	
	for (int a = 0; a < len; a++) {
	    if (i[a] > j[a]) {
		weights[i[a]][tmpdeg[i[a]]] = v[a];
		weights[j[a]][tmpdeg[j[a]]] = v[a];
		volume += 2*v[a];
	    
		// backInd[i[a]][tmpdeg[i[a]]] = tmpdeg[j[a]];
		// backInd[j[a]][tmpdeg[j[a]]] = tmpdeg[i[a]];
	    
		nbrs[i[a]][tmpdeg[i[a]]++] = j[a];
		nbrs[j[a]][tmpdeg[j[a]]++] = i[a];
	    }
	}
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
     * walk in BFS order, setting the fiel bfs
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
		    depth[nbr] = depth[curNode]+1;
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
    public int[] dfs () {

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
		dfsSub(x,c);
		c = c + 1;
	    }
	}
	
	return dfs;
    }

    private void dfsSub(int x, int c) {
	dfs[dfsPtr++] = x;
	seen[x] = true;
	comp[x] = c;
	for (int i = 0; i < deg[x]; i++) 
	    if (!seen[nbrs[x][i]]) {
		parent[nbrs[x][i]] = x;
		dfsSub(nbrs[x][i],c);
	    }

    }
	
    public int[] getDfs () {return dfs;}


    /**
     * Compute connected components using a bfs approach.
     * Returns a the characteristic vector of the components,
     * starting to index with 1
     */
    public int[] components () {

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

    
    
    /* if the graph is a tree,
       this turns it into a parent array
    */

    public int[] treeToArray() {

	dfs();
	
	int[] pArray = new int[nv];

	// set all others to parent
	for (int i = 0; i < nv; i++)
	    pArray[i] = parent[i];

	// set root to itself
	pArray[dfs[0]] = dfs[0];

	return pArray;
	
    }



    public void dump() {

	for (int x = 0 ; x < nv; x++) {
	    System.out.print(x +" : ");
	    for (int i = 0; i < deg[x]; i++) 
		System.out.print(nbrs[x][i] + " ");
	    System.out.println();
	}
	
    }


    public void compWtedDeg() {
	wtedDeg = new double[nv];
	for (int x = 0; x < nv; x++) {
	    double sum = 0;
	    for (int i = 0; i < deg[x]; i++)
		sum += weights[x][i];
	    wtedDeg[x] = sum;
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
	double[] i = new double[2*ne];
	double[] j = new double[2*ne];
	double[] v = new double[2*ne];

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
    public WeightedGraph copy () {
	WeightedGraph G = new WeightedGraph();

	G.nv = nv;
	G.ne = ne;
	G.deg = new int[nv];

	for (int a = 0; a < nv; a++) {
	    G.deg[a] = deg[a];
	}

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

    
}
