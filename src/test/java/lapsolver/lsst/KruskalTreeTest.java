package lapsolver.lsst;

import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.Stretch;
import lapsolver.generators.Grid2;
import org.junit.Test;

import static org.junit.Assert.*;

public class KruskalTreeTest {

    public static final double NS_TO_S = 1000000000.0;

    @Test
    public void testTree() throws Exception {
        final Grid2 grid = new Grid2(100);
        final SpanningTreeStrategy treeSovler = new KruskalTree();
        final Graph graph = grid.generateGraph();

        Tree lsst = null;
        final long start = System.nanoTime();
        final int nTrials = 10000;
        for (int i = 0; i < nTrials; i++) {
            lsst = treeSovler.getTree(graph);
        }
        final long elapsed = System.nanoTime() - start;
        final double time = elapsed / NS_TO_S / nTrials;

        System.out.println("Average time = " + time + " seconds.");

        final double stretch = Stretch.compute(graph, lsst).total / graph.ne;
        assertTrue("Kruskal is Low-Stretch!", stretch > 10.0);
    }
}