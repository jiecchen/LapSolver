/**
 * @file TreeSeparator.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Mon Jun 9 2014
 *
 * Find a tree separator, a vertex whose deletion splits the tree with N vertices
 * into components with at most N/2 vertices.
 *
 * Complexity: O(n)
 */

package lapsolver.algorithms;

import lapsolver.Tree;
import lapsolver.util.TreeUtils;

import java.util.Arrays;

public class TreeSeparator {
    // Find a tree separator.
    public static int find (Tree tree) {
        int[] order = TreeUtils.dfsOrder(tree);
        int[] sizes = new int[tree.nv]; // subtree sizes
        boolean[] notSep = new boolean[tree.nv]; // eliminated candidate

        // compute subtree sizes from bottom up, eliminating non-separators
        for (int i = tree.nv-1; i >= 0; i--) {
            // get subtree sizes
            int v = order[i];
            int parent = tree.getNode(v).getParent().getId();
            sizes[v] += 1;
            if (v != tree.root) {
                sizes[parent] += sizes[v];

                // too big to be child?
                if (sizes[v] > tree.nv / 2) {
                    notSep[parent] = true;
                }
            }

            // too small to be parent?
            if (tree.nv - sizes[v] > tree.nv/2) {
                notSep[v] = true;
            }
        }

        // find separator
        for (int i = 0; i < tree.nv; i++) {
            if (!notSep[i]) return i;
        }

        throw new Error("couldn't find separator");
    }
}
