package lapsolver.generators;

import lapsolver.Graph;
import lapsolver.util.NativeLoader;

public class Grid3 implements GraphFactory {
    static {
        NativeLoader.loadLibrary("lapsolver");
    }

	private final int xsize;
	private final int ysize;
	private final int zsize;
	private final double yweight;
	private final double zweight;

	private Graph graph = null;

	public Grid3(int dimension) {
		this(dimension, dimension, dimension, 1, 1);
	}

	public Grid3(int x, int y, int z) {
		this(x, y, z, 1, 1);
	}

	public Grid3(int x, int y, int z, double yw, double zw) {
		this.xsize = x;
		this.ysize = y;
		this.zsize = z;
		this.yweight = yw;
		this.zweight = zw;
	}

	public Graph generateGraph() {
		if (graph != null) {
			return graph;
		}

		// number of edges
		int ne = (xsize * ysize) * (zsize - 1) +
				 (ysize * zsize) * (xsize - 1) +
				 (zsize * xsize) * (ysize - 1);

		int[] src = new int[ne];
		int[] dst = new int[ne];
		double[] weight = new double[ne];

		populateg3(src, dst, weight, xsize, ysize, zsize, yweight, zweight);

		graph = new Graph(src, dst, weight);
		return graph;
	}

	private native void populateg3(int[] src, int[] dst, double[] weight,
                                  int x, int y, int z, double yweight, double zweight);

}


