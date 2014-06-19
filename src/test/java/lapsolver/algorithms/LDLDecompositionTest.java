package lapsolver.algorithms;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.util.EdgeListLoader;
import lapsolver.util.GraphUtils;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class LDLDecompositionTest {
    EdgeList testCaseEdges;

    @Before
    public void setUp() throws Exception {
        testCaseEdges = new EdgeListLoader().loadFromCSV("ldltest.csv");
    }

    @Test
    public void testSolve() throws Exception {
        Graph graph = new Graph(testCaseEdges);

        double[] X = new double[graph.nv];
        for (int i = 0; i < X.length; i++) {
            X[i] = Math.random();
        }

        GraphVertexRemoval.AnswerPair gvmPair = new GraphVertexRemoval(graph).solve();
        int[] perm = gvmPair.v;
        int numRemoved = gvmPair.n;

        Graph g = new Graph(GraphUtils.permuteGraph(graph, perm));
        LDLDecomposition ldl = new LDLDecomposition(g, X);
        LDLDecomposition.ReturnPair result = ldl.solve(numRemoved);
    }
}