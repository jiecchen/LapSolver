package lapsolver.generators;

import lapsolver.Graph;

public class Grid2 implements GraphFactory {
    static {
        System.loadLibrary("lapsolver");
    }

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

        //i -> from, j -> to, v -> weight
        int src[] = new int[ne];
        int dst[] = new int[ne];
        double weight[] = new double[ne];

        populate(src, dst, weight);
//        populateC(src, dst, weight, height, width, verticalWeight);

        graph = new Graph(src, dst, weight);
        return graph;
    }

    private void populate(int[] src, int[] dst, double[] weight) {
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


    private native void populateC(int[] src, int[] dst, double[] weight, int height, int width, int verticalWeight);

    // benchmark
    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        new Grid2(2000, 2000).generateGraph();
        long endTime = System.currentTimeMillis();
        System.out.println("Total execution time = " + ((endTime - startTime) / 1000.0) + "s");
    }
}
