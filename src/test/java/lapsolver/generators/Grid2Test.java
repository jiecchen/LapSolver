package lapsolver.generators;

import org.junit.Before;
import org.junit.Test;

public class Grid2Test {
    Grid2 testGrid;

    @Before
    public void setUp() throws Exception {
        testGrid = new Grid2(1000,1000);
    }

    @Test
    public void testGenerateGraph() throws Exception {
        testGrid.generateGraph();
    }
}