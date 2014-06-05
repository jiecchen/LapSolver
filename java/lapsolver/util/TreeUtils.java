/**
 * @file TreeUtils.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jun 4 2014
 *
 * Various static methods to compute useful things on trees.
 */

package lapsolver.util;

import lapsolver.Tree;

public class TreeUtils {
    // return a permutation of the vertices in BFS order
    public static int[] bfsOrder(Tree tree) {
        int[] order = new int[tree.nv];

        // start at root
        order[0] = tree.getRoot();
        int orderPtr = 1;
        int curPtr = 1;
        int curNode = tree.getRoot();

        // at each step, expand to children of curNode, then advance curNode
        while (curPtr < tree.nv) {
            for (int i = 0; i < tree.getNode(curNode).getNumberOfChildren(); i++)
                order[orderPtr++] = tree.getNode(curNode).getChild(i);
            curNode = order[curPtr++];
        }

        return order;
    }

    // return a permutation of the vertices in DFS order
    public static int[] dfsOrder(Tree tree) {
        int[] order = new int[tree.nv];
        int[] stack = new int[tree.nv];

        // start at root
        stack[0] = tree.getRoot();
        int stack_pos = 1, order_pos = 0;

        // do DFS
        while(stack_pos > 0) {
            int v = stack[--stack_pos]; // pop
            order[order_pos++] = v;

            // push children
            for(int ch : tree.getNode(v).getChildren()) {
                stack[stack_pos++] = ch;
            }
        }

        return order;
    }

    // return depths, given a valid BFS order
    public static double[] depthFromBfsOrder(Tree tree, int[] order) {
        double[] depth = new double[tree.nv];

        // at any point, we have order[i]'s parent's depth
        depth[tree.getRoot()] = 0;
        for (int i = 1; i < tree.nv; i++) {
            depth[ order[i] ] = depth[tree.getNode(order[i]).getParent().getId()]
                              + tree.getNode(order[i]).getLength();
        }

        return depth;
    }

    // return depths, given a tree (does the BFS too)
    public static double[] getDepths(Tree tree) {
        return depthFromBfsOrder( tree, bfsOrder(tree) );
    }

    // dumps the tree to sout in a BFS ordering
    public static void dumpBFSTree(Tree tree) {
        for (int i : bfsOrder(tree)) {
            System.out.println("Vertex " + i + " with parent " +
                    tree.getNode(i).getParent().getId() + " with cost " +
                    tree.getNode(i).getLength() + " on the edge to the parent");
        }
    }
}