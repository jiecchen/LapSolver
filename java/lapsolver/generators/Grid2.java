package lapsolver.generators;

import lapsolver.Graph;

public class Grid2 implements GraphFactory {
    static {
        System.loadLibrary("lapsolver");
    }

    private static boolean USE_JNI = false;

    private final int width;
    private final int height;
    private final int verticalWeight;

    private Graph graph = null;

    /**
     * Construct an evenly-weighted 2-dimensional grid graph
     *
     * @param width  width of graph in vertices
     * @param height height of graph in vertices
     */
    public Grid2(int height, int width) {
        this.width = width;
        this.height = height;
        this.verticalWeight = 1;
    }

    /**
     * Unevenly-weighted version. Weights the vertical edges according
     * to
     *
     * @param width          width of graph in vertices
     * @param height         height of graph in vertices
     * @param verticalWeight adjusted weight for vertical edges
     */
    public Grid2(int height, int width, int verticalWeight) {
        this.width = width;
        this.height = height;
        this.verticalWeight = verticalWeight;
    }

    private int getIdx(final int i, final int j) {
        return width * i + j;
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

        int [] srcArr = new int[ne];
        int [] dstArr = new int[ne];
        double [] weightArr = new double[ne];

        if(USE_JNI) {
            populateC(srcArr, dstArr, weightArr, height, width, verticalWeight);
        } else {
            populate(srcArr, dstArr, weightArr);
        }

        graph = new Graph(srcArr, dstArr, weightArr);
        return graph;
    }

    /**
     * Native C version of populateC below. Not much faster than the JIT version.
     * @param src
     * @param dst
     * @param weight
     * @param height
     * @param width
     * @param verticalWeight
     */
    private native void populateC(int[] src, int[] dst, double[] weight,
                                  int height, int width, int verticalWeight);

    /**
     * Populate the src, dst, and weight arrays with the grid edges
     *
     * @param src    the source vertices
     * @param dst    the corresponding destinations
     * @param weight the weight of each edge
     */
    private void populate(int[] src, int[] dst, double[] weight) {
        // populate edge lists
        int e = 0;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int cur = getIdx(i, j);
                if (i + 1 > 0 && i + 1 < height) { //vertical edge
                    src[e] = cur;
                    dst[e] = getIdx(i + 1, j);
                    weight[e++] = verticalWeight;
                }
                if (j + 1 > 0 && j + 1 < width) { //horizontal edge
                    src[e] = cur;
                    dst[e] = getIdx(i, j + 1);
                    weight[e++] = 1;
                }
            }
        }
    }

    // benchmark
    public static void main(String[] args) {
        int nBenchmarks = 10;
        int size = 500;

        for (int i = 0; i < 10; i++) {
            System.out.print("Without JNI, ");
            runBenchmark(nBenchmarks, false, size);

            System.out.print("With JNI, ");
            runBenchmark(nBenchmarks, true, size);

            System.out.println("");
        }
    }

    private static void runBenchmark(int nBenchmarks, boolean jni, int size) {
        USE_JNI = jni;
        long startTime = System.currentTimeMillis();
        for (int i = 0; i < nBenchmarks; i++)
            new Grid2(size, size).generateGraph();
        long endTime = System.currentTimeMillis();
        System.out.println("average execution time = " + ((endTime - startTime) / (1000.0 * nBenchmarks)) + "s");
    }
}
