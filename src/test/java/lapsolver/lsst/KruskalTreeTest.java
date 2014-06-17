package lapsolver.lsst;

import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.Stretch;
import lapsolver.generators.Grid2;
import org.junit.Test;

import static org.junit.Assert.*;

public class KruskalTreeTest {

    @Test
    public void testTree() throws Exception {
        final Grid2 grid = new Grid2(100);
        final SpanningTreeStrategy treeSovler = new KruskalTree();
        final Graph graph = grid.generateGraph();
        final Tree lsst = treeSovler.getTree(graph);
        final double stretch = Stretch.compute(graph, lsst).total / graph.ne;

        assertTrue("Kruskal is Low-Stretch!", stretch > 10.0);
    }
}