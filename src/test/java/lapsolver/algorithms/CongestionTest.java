/**
 * @file CongestionTest.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 14, 2014
 *
 * A simple test of the Congestion algorithm for a weighted tree.
 */

package lapsolver.algorithms;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

public class CongestionTest {
    @Test
    public void testCompute() throws Exception {
        // build tree
        int[] binTree_p = {0,0,0,1,1,2,2};
        double[] binTree_w = {-1,4,4,1,1,1,1};
        Tree binTree = new Tree(binTree_p, binTree_w);
        Assert.assertEquals(0, binTree.root);

        // build clique
        EdgeList clique_edges = new EdgeList(21);
        int edgePos = 0;
        for (int i = 0; i < 7; i++) {
            for (int j = 0; j < i; j++) {
                clique_edges.u[edgePos] = i;
                clique_edges.v[edgePos] = j;
                clique_edges.weight[edgePos] = 1;
                edgePos++;
            }
        }
        Graph clique = new Graph(clique_edges);

        // get and check congestion tree
        Tree congTree = Congestion.compute(clique, binTree);
        EdgeList edges = new EdgeList(congTree);
        double[] expect = {3.0, 3.0, 6.0, 6.0, 6.0, 6.0};
        Assert.assertArrayEquals(expect, edges.weight, 1e-8);
    }
}