package lapsolver.generators;

import lapsolver.Graph;
import lapsolver.util.NativeLoader;

public class Hypercube implements GraphFactory {
	static {
		NativeLoader.loadLibrary("lapsolver");
	}

	private final int wsize;
	private final int xsize;
	private final int ysize;
	private final int zsize;
	private final double xweight;
	private final double yweight;
	private final double zweight;

	private Graph graph = null;

	public Hypercube(int dimension) {
		this(dimension, dimension, dimension, dimension, 1.0, 1.0, 1.0);
	}

	public Hypercube(int w, int x, int y, int z) {
		this(w, x, y, z, 1.0, 1.0, 1.0);
	}

	public Hypercube(int w, int x, int y, int z, double xw, double yw, double zw) {
		this.wsize = w;
		this.xsize = x;
		this.ysize = y;
		this.zsize = z;
		this.xweight = xw;
		this.yweight = yw;
		this.zweight = zw;
	}

	public Graph generateGraph() {
		if (graph != null) {
			return graph;
		}

		// number of edges
		int ne = (xsize * ysize * zsize) * (wsize - 1) +
				 (wsize * ysize * zsize) * (xsize - 1) +
				 (wsize * xsize * zsize) * (ysize - 1) +
				 (wsize * xsize * ysize) * (zsize - 1);

		int[] src = new int[ne];
		int[] dst = new int[ne];
		double[] weight = new double[ne];

		populateg4(src, dst, weight, wsize, xsize, ysize, zsize,
				   xweight, yweight, zweight);

		graph = new Graph(src, dst, weight);
		return graph;
	}

	private native void populateg4(int[] src, int[] dst, double[] weight,
								   int w, int x, int y, int z, double xweight,
								   double yweight, double zweight);

}