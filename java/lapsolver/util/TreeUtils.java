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
        order[0] = tree.root;
        int orderPtr = 1;
        int curPtr = 1;
        int curNode = tree.root;

        // at each step, expand to children of curNode, then advance curNode
        while (curPtr < tree.nv) {
            for (int i = 0; i < tree.nodes[curNode].getNumberOfChildren(); i++)
                order[orderPtr++] = tree.nodes[curNode].getChild(i);
            curNode = order[curPtr++];
        }

        return order;
    }

    // return depths, given a valid BFS order
    public static double[] depthFromBfsOrder(Tree tree, int[] order) {
        double[] depth = new double[tree.nv];

        // at any point, we have order[i]'s parent's depth
        depth[tree.root] = 0;
        for (int i = 1; i < tree.nv; i++) {
            depth[ order[i] ] = depth[ tree.nodes[order[i]].parent ] + tree.nodes[order[i]].length;
        }

        return depth;
    }
}
