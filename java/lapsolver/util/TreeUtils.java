/**
 * @file TreeUtils.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jun 4 2014
 *
 * Various static methods to compute useful things on trees.
 */

package lapsolver.util;

import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.EdgeList;

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

    // return depths, given a valid BFS or DFS tree order
    public static double[] depthFromTreeOrder(Tree tree, int[] order) {
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
        return depthFromTreeOrder(tree, bfsOrder(tree));
    }

    // dumps the tree to sout in a BFS ordering
    public static void dumpBFSTree(Tree tree) {
        for (int i : bfsOrder(tree)) {
            System.out.println("Vertex " + i + " with parent " +
                    tree.getNode(i).getParent().getId() + " with cost " +
                    tree.getNode(i).getLength() + " on the edge to the parent");
        }
    }

    // EdgeList from a spanning tree: return off-tree edges
    public EdgeList offTreeEdges (Graph graph, Tree spanningTree) {
        int ne = graph.ne - graph.nv + 1;

        int[] u = new int[ne];
        int[] v = new int[ne];
        double[] w = new double[ne];

        int index = 0;
        int vertexCount = graph.nv;
        for (int i = 0; i < vertexCount; i++) {
            for (int j = 0; j < graph.deg[i]; j++) {
                int dest = graph.nbrs[i][j];
                double wt = graph.weights[i][j];

                // skip tree edges
                if (spanningTree.getNode(i).getParent().getId() == dest ||
                        spanningTree.getNode(dest).getParent().getId() == i) {
                    continue;
                }

                // only count an edge once
                if (i < dest) {
                    u[index] = i;
                    v[index] = dest;
                    w[index] = wt;
                    index++;
                }
            }
        }

        return new EdgeList(u, v, w);
    }

}
