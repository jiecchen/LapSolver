package lapsolver.generators;

import lapsolver.Graph;
import lapsolver.util.NativeLoader;

public class CompleteBinaryTree implements GraphFactory {
	static {
		NativeLoader.loadLibrary("lapsolver");
	}

	private final int h;
	private Graph graph = null;

	public CompleteBinaryTree(int height) {
		this.h = height;
	}

	@Override
	public Graph generateGraph() {
		if (graph != null) {
			return graph;
		}

		int ne = (2 << h) - 2;

		int[] src = new int[ne];
		int[] dst = new int[ne];
		double[] weight = new double[ne];

		populateCBT(src, dst, weight, h);

		graph = new Graph(src, dst, weight);
		return graph;
	}

	private native void populateCBT(int[] src, int[] dst,
								    double[] weight, int h);
}