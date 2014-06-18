package lapsolver.generators;

import lapsolver.Graph;
import lapsolver.util.NativeLoader;

public class TriangleGrid2 implements GraphFactory {
    static {
        NativeLoader.loadLibrary("lapsolver");
    }

    private final int width;
    private final int height;
    private final int verticalWeight;

    private Graph graph = null;

    /**
     * Construct a square, evenly-weighted 2-dimensional grid graph
     *
     * @param dimension both height and width in vertices
     */
    public TriangleGrid2(int dimension) {
        this(dimension, dimension);
    }

    /**
     * Construct an evenly-weighted 2-dimensional grid graph
     *
     * @param width  width of graph in vertices
     * @param height height of graph in vertices
     */
    public TriangleGrid2(int height, int width) {
        this(height, width, 1);
    }

    /**
     * Unevenly-weighted version. Weights the vertical edges according
     * to
     *
     * @param width          width of graph in vertices
     * @param height         height of graph in vertices
     * @param verticalWeight adjusted weight for vertical edges
     */
    public TriangleGrid2(int height, int width, int verticalWeight) {
        this.width = width;
        this.height = height;
        this.verticalWeight = verticalWeight;
    }

    /**
     * Actually generate the graph
     *
     * @return a Graph representing a 2D grid
     */
    @Override
    public Graph generateGraph() {
        if (graph != null)
            return graph;

        //number of edges, vertices, non-bottom row
        int ne = (2 * width * height) - width - height;
        ne += (width - 1) * (height - 1);

        int[] src = new int[ne];
        int[] dst = new int[ne];
        double[] weight = new double[ne];

        populateC(src, dst, weight, height, width, verticalWeight);

        graph = new Graph(src, dst, weight);
        return graph;
    }

    /**
     * Native C call to vectorize the grid construction
     *
     * @param src            the source vertices
     * @param dst            the corresponding destinations
     * @param weight         the weight of each edge
     * @param height         same as the corresponding field in the class
     * @param width          same as the corresponding field in the class
     * @param verticalWeight same as the corresponding field in the class
     */
    private native void populateC(int[] src, int[] dst, double[] weight,
                                  int height, int width, int verticalWeight);

}
