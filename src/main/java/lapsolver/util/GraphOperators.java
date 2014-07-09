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
    public static RealLinearOperator buildLaplacianOperator(final Graph graph) {
        return buildLaplacianOperator(graph, new double[graph.nv]);
    }

    public static RealLinearOperator buildLaplacianOperator(final Graph graph, final double[] d) {
        return new RealLinearOperator() {
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
                if (graph.nv != x.getDimension())
                    throw new DimensionMismatchException(x.getDimension(), graph.nv);
                double[] xarr = x.toArray();
                double[] lx = GraphUtils.applyLaPlacian(graph, xarr);
                for (int i = 0; i < graph.nv; i++)
                    lx[i] += d[i] * xarr[i];
                return new ArrayRealVector(lx);
            }
        };
    }
}
