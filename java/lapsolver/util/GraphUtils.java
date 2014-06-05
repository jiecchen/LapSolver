/**
 * @file GraphUtils.java
 * @author Alex Reinking <alexander.reinking@gmail.com>
 * @date Thu Jun 5 2014
 *
 * Various static methods to compute useful things on graphs.
 */

package lapsolver.util;

import lapsolver.Graph;
import lapsolver.Tree;

public class GraphUtils {
    // given a graph structure that represents a tree, turn it into a tree
    public Tree toTree(Graph g) {
        return g.treeToTree();
    }

    // prints the contents of g to stdout, as an adjacency list
    public void dump(Graph g) {
        for (int x = 0; x < g.nv; x++) {
            System.out.print(x + " : ");
            for (int i = 0; i < g.deg[x]; i++)
                System.out.print(g.nbrs[x][i] + " ");
            System.out.println();
        }
    }
}
