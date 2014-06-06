/**
 * @file Congestion.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * Dan's old code for computing stretch.
 */

package lapsolver.algorithms;

import lapsolver.Tree;
import lapsolver.Graph;
import lapsolver.util.TreeUtils;

import java.util.ArrayDeque;

/**
 * Created by cyrilzhang on 6/5/14.
 */
public class StretchDan {
    Graph graph, stretchGraph;
    Tree spanningTree;

    ArrayDeque<Integer> aux;
    ArrayDeque<Integer> aux2;
    int[] comp;
    int[] order;
    double[] depth;

    public StretchDan(Graph graph, Tree spanningTree) {
        this.graph = graph;
        this.spanningTree = spanningTree;
    }

    // compute total stretch of spanningTree in graph
    public double compTotalStretch() {
        aux = new ArrayDeque<>(graph.nv);
        aux2 = new ArrayDeque<>(graph.nv);

        order = TreeUtils.bfsOrder(spanningTree);
        depth = TreeUtils.depthFromTreeOrder(spanningTree, order);
        comp = new int[graph.nv]; // the component of a node

        // the algorithm will create a new component
        // each time it encounters a node with no children

        int[] nodeSize = new int[graph.nv];  // how much is below a node

        int orderPtr; // pointer into the order
        int node;     // and, its node
        int numComps = 0;

        double totStretch = 0;

        for (orderPtr = graph.nv-1; orderPtr >= 0; orderPtr--) {
            node = order[orderPtr];

            if (spanningTree.getNode(node).getNumberOfChildren() == 0) {
                comp[node] = ++numComps;
                nodeSize[node] = graph.deg[node];
                // nothing to merge
            } else {
                // find largest kid, and set this node size
                nodeSize[node] = graph.deg[node];

                int maxKid = -1;
                int maxKidSize = -1;

                for (int kid : spanningTree.getNode(node).getChildren()) {
                    int kidSize = nodeSize[kid];

                    nodeSize[node] += kidSize;

                    if (kidSize > maxKidSize) {
                        maxKidSize = kidSize;
                        maxKid = kid;
                    }
                }

                int curComp = comp[maxKid];

                // go through other children,
                // re-componentizing,
                // and checking edges for new merge (if other end is this comp)

                for (int kid : spanningTree.getNode(node).getChildren()) {
                    if (kid != maxKid)
                        totStretch += totStretchTraverse(kid, node, curComp);
                }

                // and, still have to fix up the root!
                // note that this code is easy,
                // but costs a branch in totStretchTraverse
                // that is unnecessary
                // would be faster to handle this node
                // on its own (because do not go beneath)
                totStretch += totStretchTraverse(node, node, curComp);
            }
        }
        return totStretch;
    }

    /*
     * the subroutine for compStretch
     *
     * first, account for all edges between this subtree
     *  and those in comp curComp
     *
     * then, make the component of all these curComp
     */
    public double totStretchTraverse(int child, int node, int curComp) {
        double stretch = 0;

        // now, traverse
        aux.clear();
        aux.add(child);

        while (!aux.isEmpty()) {
            int curNode = aux.poll();

            // note: following test is unnecessary
            // except when we are handling the root of the tree

            if (comp[curNode] != curComp) {
                for (int v : spanningTree.getNode(node).getChildren())
                    aux.add(v);

                // put node on list to change comp
                aux2.add(curNode);

                // now, go over the edges to check if can now compute stretch
                for (int i = 0; i < graph.deg[curNode]; i++) {
                    if (comp[graph.nbrs[curNode][i]] == curComp) {
                        int nbr = graph.nbrs[curNode][i];
                        double len = 1 / graph.weights[curNode][i];
                        stretch += (depth[curNode] + depth[nbr] - 2 * depth[node]) / len;
                        /*
                        System.out.println("(" + curNode + ", " + nbr + ") : " +
                                   + (depth[curNode] + depth[nbr] - 2*depth[node]));
                        */
                    }
                }
            }
        }

        //--------------------------------------------
        //
        // finally, recomp all the nodes we've seen

        while (!aux2.isEmpty()) {
            int curNode = aux2.poll();
            comp[curNode] = curComp;
        }

        return stretch;
    }

    /**
     * Computes the stretch of every edge
     *
     * @return a weighted graph where the weight is the stretch of every edge
     */
    public Graph compStretches() {
        stretchGraph = graph.copy();

        // make all stretches zero, just to find errors
        for (int a = 0; a < graph.nv; a++) {
            for (int i = 0; i < graph.deg[a]; i++) {
                stretchGraph.weights[a][i] = 0;
            }
        }

        order = TreeUtils.bfsOrder(spanningTree);
        depth = TreeUtils.depthFromTreeOrder(spanningTree, order);
        comp = new int[graph.nv];   // the component of a node
        //
        // the algorithm will create a new component
        // each time it encounters a node with no children

        int[] nodeSize = new int[graph.nv];  // how much is below a node

        int orderPtr; // pointer into the order
        int node;     // and, its node
        int numComps = 0;

        aux = new ArrayDeque<>(graph.nv);  // will use for traversing
        aux2 = new ArrayDeque<>(graph.nv);  // will use for traversing

        double totStretch = 0;

        for (orderPtr = graph.nv-1; orderPtr >= 0; orderPtr--) {
            node = order[orderPtr];

            if (spanningTree.getNode(node).getNumberOfChildren() == 0) {
                comp[node] = ++numComps;
                nodeSize[node] = graph.deg[node];
                // nothing to merge
            } else {

                // find largest kid, and set this node size
                nodeSize[node] = graph.deg[node];

                int maxKid = -1;
                int maxKidSize = -1;

                for (int kid : spanningTree.getNode(node).getChildren()) {
                    int kidSize = nodeSize[kid];

                    nodeSize[node] += kidSize;

                    if (kidSize > maxKidSize) {
                        maxKidSize = kidSize;
                        maxKid = kid;
                    }
                }

                int curComp = comp[maxKid];

                // go through other children,
                // re-componentizing,
                // and checking edges for new merge (if other end is this comp)

                for (int kid : spanningTree.getNode(node).getChildren()) {
                    if (kid != maxKid)
                        totStretch += allStretchTraverse(kid, node, curComp);
                }

                // and, still have to fix up the root!
                // note that this code is easy,
                // but costs a branch in allStretchTraverse
                // that is unnecessary
                // would be faster to handle this node
                // on its own (because do not go beneath)
                totStretch += allStretchTraverse(node, node, curComp);
            }
        }

        return stretchGraph;
    }

    /*
     * the subroutine for compStretch
     *
     * first, account for all edges betwen this subtree
     *  and those in comp curComp
     *
     * then, make the component of all these curComp
     */
    public double allStretchTraverse(int kid, int node, int curComp) {
        double stretch = 0;

        // now, traverse
        aux.clear();
        aux.add(kid);

        while (!aux.isEmpty()) {
            int curNode = aux.poll();

            // note: following test is unnecessary
            // except when we are handling the root of the tree

            if (comp[curNode] != curComp) {
                for (int v : spanningTree.getNode(node).getChildren())
                    aux.add(v);

                // put node on list to change comp
                aux2.add(curNode);

                // now, go over the edges to check if can now compute stretch
                for (int i = 0; i < graph.deg[curNode]; i++) {
                    if (comp[graph.nbrs[curNode][i]] == curComp) {
                        int nbr = graph.nbrs[curNode][i];
                        double len = 1 / graph.weights[curNode][i];
                        double thisStretch = (depth[curNode] + depth[nbr] - 2 * depth[node]) / len;
                        stretch += thisStretch;
                        stretchGraph.weights[curNode][i] = thisStretch;

                /*
                System.out.println("(" + curNode + ", " + nbr + ") : " +
                           + (depth[curNode] + depth[nbr] - 2*depth[node]));
                */
                    }
                }
            }
        }

        //--------------------------------------------
        //
        // finally, recomp all the nodes we've seen
        while (!aux2.isEmpty()) {
            int curNode = aux2.poll();
            comp[curNode] = curComp;
        }

        return stretch;
    }
}
