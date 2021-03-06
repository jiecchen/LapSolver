/**
 * @file Tree.java
 * @author Daniel Spielman <spielman@cs.yale.edu>
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jan 8 2014
 *
 * (c) Daniel Spielman, 2014, part of the YINSmat package
 *
 * A static tree data structure.
 * Heavily refactored by Cyril on Thursday, June 12, 2014.
 *
 */

package lapsolver;

import lapsolver.util.GraphUtils;
import lapsolver.util.TreeUtils;

public class Tree {
    public int nv;    // number of vertices (ne = nv-1)
    public int root;  // auxiliarySize of root of tree
    public int[] parent;     // parent[permutation] = id of permutation's parent, parent[root] = root
    public double[] weight;  // weight[permutation] = edge weight of (permutation, parent[permutation])
    public int[][] children; // children[permutation][0..nChildren[permutation]-1] = indices of children

    // tree from parent array and weights
    public Tree(int[] parent, double[] weight) {
        nv = parent.length;
        this.parent = parent;
        this.weight = weight;
        int[] nChildren = new int[nv];

        if (weight == null) {
            this.weight = new double[nv];
            for (int i = 0; i < nv; i++) {
                this.weight[i] = 1.0;
            }
        }

        for (int i = 0; i < nv; i++) {
            if (parent[i] == i) { // find root
                root = i;
            }
            else { // count child
                nChildren[parent[i]]++;
            }
        }

        // build children arrays from parents
        children = new int[nv][];
        for (int i = 0; i < nv; i++) {
            children[i] = new int[nChildren[i]];
        }

        int[] childPos = new int[nv];
        for (int i = 0; i < nv; i++) {
            if (i != root) {
                children[parent[i]][childPos[parent[i]]++] = i;
            }
        }
    }

    // build tree from parent array, setting all lengths to 1
    public Tree(int[] parent) {
        this(parent, null);
    }

    // build tree from undirected acyclic graph
    public Tree(Graph treeGraph) {
        this(GraphUtils.toTree(treeGraph));
    }

    public Tree(EdgeList edges) { this(new Graph(edges)); }

    // copy constructor
    public Tree (Tree other) {
        nv = other.nv;
        root = other.root;
        parent = other.parent.clone();
        weight = other.weight.clone();
        children = new int[nv][];

        for (int i = 0; i < nv; i++) {
            children[i] = other.children[i].clone();
        }
    }
}
