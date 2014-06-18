package lapsolver.generators;

import lapsolver.Graph;

import java.util.Arrays;

public class PetersenGraph implements GraphFactory {
    private int n;
    private int k;

    /**
     * Create a generalized Petersen Graph with parameters n and k.
     *
     * The classic Petersen Graph has n=5, k=2
     *
     * @param n The number of points on the inner star
     * @param k Determines the shape of the inner star. 1 <= k < n/2
     */
    public PetersenGraph(int n, int k) {
        if (k < 1 || k >= Math.ceil((double) n / 2))
            throw new IllegalArgumentException("k must be in [1, n/2)");

        this.n = n;
        this.k = k;
    }

    @Override
    public Graph generateGraph() {
        int[] src = new int[3 * n];
        int[] dst = new int[3 * n];
        double[] weights = new double[3 * n];

        Arrays.fill(weights, 1.0); // TODO: customize the weights

        int e = 0;
        for (int i = 0; i < n; i++) {
            src[e] = i;
            dst[e++] = (i + 1) % n;

            src[e] = i;
            dst[e++] = n + i;

            src[e] = n + i;
            dst[e++] = n + ((i + k) % n);
        }

        return new Graph(src, dst, weights);
    }
}
