/**
 * @file GraphOperators.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jul 2 2014
 *
 * Wrappers for Apache Commons linear operators.
 */

package lapsolver.util;

import lapsolver.Graph;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.apache.commons.math3.linear.RealVector;

public class GraphOperators {
    public static class LaplacianOperator extends RealLinearOperator {
        public Graph graph;
        public double[] d;

        public LaplacianOperator (Graph graph) {
            this.graph = graph;
            d = new double[graph.nv];
        }

        public LaplacianOperator (Graph graph, double[] d) {
            this.graph = graph;
            this.d = d;
        }

        @Override
        public int getRowDimension() {
            return graph.nv;
        }

        @Override
        public int getColumnDimension() {
            return graph.nv;
        }

        @Override
        public RealVector operate(RealVector x) throws DimensionMismatchException {
            double[] xarr = x.toArray();
            double[] lx = GraphUtils.applyLaplacian(graph, xarr);
            for (int i = 0; i < graph.nv; i++) {
                lx[i] += d[i] * xarr[i];
            }
            return new ArrayRealVector( lx );
        }
    }
}
