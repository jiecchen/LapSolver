package lapsolver.generators;

import lapsolver.WeightedGraph;

import java.util.Arrays;

public class Grid2 implements GraphFactory {
    private final int width;
    private final int height;
    private final int verticalWeight;

    private WeightedGraph graph = null;

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

    private int getIdx(int i, int j) {
        return width * i + j;
    }

    /**
     * Actually generate the graph
     *
     * @return a WeightedGraph representing a 2D grid
     */
    @Override
    public WeightedGraph generateGraph() {
        if (graph != null)
            return graph;

        //number of edges, vertices, non-bottom row
        int ne = (2 * width * height) - width - height;

        //i -> from, j -> to, v -> weight
        int src[] = new int[ne];
        int dst[] = new int[ne];
        double weight[] = new double[ne];
        Arrays.fill(weight, 1);

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

        graph = new WeightedGraph(src, dst, weight);
        return graph;
    }

    // benchmark
//    public static void main(String[] args) {
//        System.out.print("Before: ");
//        System.out.println((double) Runtime.getRuntime().totalMemory() / 1024.0 / 1024.0);
//
//        GraphFactory gf = new Grid2(400,400);
//        double[][] ijv = gf.generateGraph().toIJV();
//
//        System.out.print("After: ");
//        System.out.println((double) Runtime.getRuntime().totalMemory() / 1024.0 / 1024.0);
//        System.out.println("Done (" + ijv[1][1] + ")");
//    }
}
