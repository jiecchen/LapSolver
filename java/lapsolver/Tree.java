/**
 * @file Tree.java
 * @author Daniel Spielman <spielman@cs.yale.edu>
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jan 8 2014
 *
 * (c) Daniel Spielman, 2014, part of the YINSmat package
 *
 * A static tree data structure.
 */

package lapsolver;

import java.util.ArrayList;

/**
 * a class enabling tree operations
 * that are useful for preconditioning
 * Fat because it is implemented for easy code,
 * not speed
 */
public class Tree {
    public int nv;
    private TreeNode[] nodes;
    private int root;

    /**
     * Make an empty Tree
     *
     * @param nv number of nodes
     */
    public Tree(int nv) {
        nodes = new TreeNode[nv];
    }

    /**
     * Make a Tree from parent array
     *
     * @param p the parent array
     */
    public Tree(int[] p) {
        nodes = new TreeNode[p.length];
        this.fromParentArray(p);
    }

    public TreeNode getNode(int n) {
        return nodes[n];
    }

    /**
     * Builds tree from parent pointers.
     * the root should have a self-loop.
     *
     * @param p pointers to parents
     */
    private void fromParentArray(int[] p) {
        nv = p.length;

        nodes = new TreeNode[nv];

        for (int i = 0; i < nv; i++) {
            nodes[i] = new TreeNode(p[i], i);
            if (p[i] == i)
                root = i;
        }

        buildChildrenFromParents();
    }

    /**
     * Assuming just parents are set,
     * fills in all the children.
     */
    public void buildChildrenFromParents() {
        for (int i = 0; i < nv; i++)
            if (i != root)
                nodes[nodes[i].parent].addChild(i);
    }

    /**
     * Make a Tree from parent array and weights on edges
     *
     * @param p the parent array
     * @param w the weight array
     */
    public Tree(int[] p, double[] w) {
        nodes = new TreeNode[p.length];
        this.fromParentLengthArray(p, w);
    }

    /**
     * Builds tree from parent pointers.
     * the root should have a self-loop.
     *
     * @param p pointers to parents
     * @param w weights on edges (which now turn to 1/length)
     */
    public void fromParentLengthArray(int[] p, double[] w) {
        nv = p.length;

        nodes = new TreeNode[nv];

        for (int i = 0; i < nv; i++) {
            nodes[i] = new TreeNode(p[i], i, 1 / w[i]);
            if (p[i] == i)
                root = i;
        }

        buildChildrenFromParents();
    }

    public int getRoot() {
        return root;
    }

    /**
     * Outputs a parent pointer array for the tree
     * The root points to itself
     *
     * @return a parent pointer array
     */
    public int[] toParentArray() {
        int[] p = new int[nv];

        for (int i = 0; i < nv; i++) {
            p[i] = nodes[i].parent;
        }

        return p;
    }

    /*
     * A node of the tree
     * default length is 1
     * length is length of edge to the parent
     */
    public class TreeNode {
        private int parent;
        private final int id;
        private double length;

        public ArrayList<Integer> getChildren() {
            return children;
        }

        private final ArrayList<Integer> children;

        /**
         * Create a node by specifying its parent
         *
         * @param p parent
         */
        public TreeNode(int p, int self) {
            parent = p;
            id = self;
            length = 1;
            children = new ArrayList<Integer>();
        }

        public TreeNode(int p, int self, double edgeLength) {
            parent = p;
            id = self;
            length = edgeLength;
            children = new ArrayList<Integer>();
        }

        public TreeNode getParent() {
            return nodes[parent];
        }

        public void setParent(TreeNode p) {
            parent = p.getId();
        }

        public double getLength() {
            return length;
        }

        public void setLength(double length) {
            this.length = length;
        }

        public void addChild(int k) {
            children.add(k);
        }

        public int getNumberOfChildren() {
            return children.size();
        }

        public int getChild(int i) {
            return children.get(i);
        }

        public int getId() {
            return id;
        }
    }
}
