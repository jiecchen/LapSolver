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

    // empty tree
    public Tree(int nv) {
        nodes = new TreeNode[nv];
    }

    // tree from parent array
    public Tree(int[] p) {
        nv = p.length;
        nodes = new TreeNode[nv];

        for (int i = 0; i < nv; i++) {
            nodes[i] = new TreeNode(p[i], i);
            if (p[i] == i)
                root = i;
        }

        buildChildrenFromParents();
    }

    // copy constructor
    public Tree (Tree other) {
        nv = other.nv;
        root = other.root;
        nodes = new TreeNode[nv];

        for (int i = 0; i < nv; i++) {
            nodes[i] = new TreeNode(other.nodes[i]);
        }
    }

    public TreeNode getNode(int n) {
        return nodes[n];
    }
    public int getRoot() {
        return root;
    }

    // assuming parents are set, build lists of children for each node
    public void buildChildrenFromParents() {
        for (int i = 0; i < nv; i++)
            if (i != root)
                nodes[nodes[i].parent].addChild(i);
    }

    // build tree from parent array with weights
    // sets lengths to reciprocals of weights
    public Tree(int[] p, double[] w) {
        nodes = new TreeNode[p.length];
        nv = p.length;
        nodes = new TreeNode[nv];

        for (int i = 0; i < nv; i++) {
            nodes[i] = new TreeNode(p[i], i, 1 / w[i]);
            if (p[i] == i)
                root = i;
        }

        buildChildrenFromParents();
    }

    // a vertex in the tree, which contain parent and children pointers
    // length is reciprocal of weight of edge from node to parent
    public class TreeNode {
        private int parent;
        private final int id;
        private double length;
        private final ArrayList<Integer> children;

        // constructor with weight 1
        public TreeNode(int parent, int id) {
            this(parent, id, 1.0);
        }

        // parent id, self id, edge length (like a directed EdgeList element)
        // remember length = 1/weight
        public TreeNode(int parent, int id, double length) {
            this.parent = parent;
            this.id = id;
            this.length = length;
            children = new ArrayList<>();
        }

        // copy constructor
        public TreeNode(TreeNode other) {
            parent = other.parent;
            id = other.id;
            length = other.length;
            children = new ArrayList<>(other.children);
        }

        public ArrayList<Integer> getChildren() {
            return children;
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
