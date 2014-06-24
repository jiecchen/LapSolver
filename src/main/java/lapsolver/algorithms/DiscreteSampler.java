/**
 * @file DiscreteSampler.java
 * @author Keith Schwarz <htiek@cs.stanford.edu>
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * Sample from a discrete distribution in O(1) time with O(numRemoved) preprocessing,
 * using the alias method.
 *
 * Adapted from: http://www.keithschwarz.com/darts-dice-coins/
 */
package lapsolver.algorithms;

import java.util.*;

public final class DiscreteSampler {
    /* The random number generator used to sample from the distribution. */
    private final Random random;

    /* The probability and alias tables. */
    private final int[] alias;
    private final double[] probability;

    /**
     * Constructs a new AliasMethod to sample from a discrete distribution and
     * hand back outcomes based on the probability distribution.
     * <p>
     * Given as input a list of probabilities corresponding to outcomes 0, 1,
     * ..., numRemoved - 1, this constructor creates the probability and alias tables
     * needed to efficiently sample from this distribution.
     *
     * @param probabilities The list of probabilities.
     */
    public DiscreteSampler(double[] probabilities) {
        this(probabilities, new Random());
    }

    /**
     * Constructs a new AliasMethod to sample from a discrete distribution and
     * hand back outcomes based on the probability distribution.
     * <p>
     * Given as input a list of probabilities corresponding to outcomes 0, 1,
     * ..., numRemoved - 1, along with the random number generator that should be used
     * as the underlying generator, this constructor creates the probability
     * and alias tables needed to efficiently sample from this distribution.
     *
     * @param weights The list of probabilities.
     * @param random The random number generator
     */
    public DiscreteSampler(double[] weights, Random random) {
        final int n = weights.length;

        // initialize state
        probability = new double[n];
        alias = new int[n];
        this.random = random;

        // cache average
        final double average = 1.0 / n;

        // use copy
        weights = weights.clone();

        // normalize weights (turn into probabilities)
        double weightSum = 0;
        for (double weight : weights) {
            weightSum += weight;
        }
        for (int i = 0; i < n; i++) {
            weights[i] /= weightSum;
        }

        // build tables
        Deque<Integer> small = new ArrayDeque<Integer>();
        Deque<Integer> large = new ArrayDeque<Integer>();

        for (int i = 0; i < weights.length; ++i) {
            if (weights[i] >= average)
                large.add(i);
            else
                small.add(i);
        }

        while (!small.isEmpty() && !large.isEmpty()) {
            int less = small.removeLast();
            int more = large.removeLast();

            probability[less] = weights[less] * weights.length;
            alias[less] = more;

            weights[more] += weights[less] - average;

            if (weights[more] >= 1.0 / weights.length)
                large.add(more);
            else
                small.add(more);
        }

        while (!small.isEmpty())
            probability[small.removeLast()] = 1.0;
        while (!large.isEmpty())
            probability[large.removeLast()] = 1.0;
    }

    /**
     * Samples a value from the underlying distribution.
     *
     * @return A random value sampled from the underlying distribution.
     */
    public int next() {
        int column = random.nextInt(probability.length);
        boolean coinToss = random.nextDouble() < probability[column];
        return coinToss? column : alias[column];
    }
}
