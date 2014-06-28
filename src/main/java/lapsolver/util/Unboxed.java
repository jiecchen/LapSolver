package lapsolver.util;

import java.util.Collection;

/**
 * Created by alexreinking on 6/27/14.
 */
public class Unboxed {
    public static int[] intsToArray(Collection<Integer> ints) {
        Integer[] intArr = new Integer[ints.size()];
        ints.toArray(intArr);
        return toPrimitive(intArr);
    }

    public static double[] doublesToArray(Collection<Double> doubles) {
        Double[] doubleArr = new Double[doubles.size()];
        doubles.toArray(doubleArr);
        return toPrimitive(doubleArr);
    }

    public static double[] toPrimitive(Double[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return new double[0];
        }
        final double[] result = new double[array.length];
        for (int i = 0; i < array.length; i++)
            result[i] = array[i];
        return result;
    }

    public static int[] toPrimitive(Integer[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return new int[0];
        }
        final int[] result = new int[array.length];
        for (int i = 0; i < array.length; i++)
            result[i] = array[i];
        return result;
    }
}
