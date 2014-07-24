/**
 * @file HighStretchSampler.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jul 24 2014
 *
 * Produces a distribution that gives probability 1 to the q edges of
 * highest stretch, and 0 to all others.
 */

package lapsolver.algorithms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class HighStretchSampler {
    public static double[] compute (final double[] allStretches, int q) {
        ArrayList<Integer> indices = new ArrayList<>(allStretches.length);
        for (int i = 0; i < allStretches.length; i++) {
            indices.add(i);
        }

        Collections.sort(indices, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return Double.compare(allStretches[o2], allStretches[o1]);
            }
        });

        double[] p = new double[allStretches.length];
        for (int i = 0; i < q; i++) {
            p[indices.get(i)] = 1;
        }

        return p;
    }
}
