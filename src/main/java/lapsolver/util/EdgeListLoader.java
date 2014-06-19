/**
 * @file StarDecompositionTree.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Thu Jun 19 2014
 *
 * Loads and constructs edge lists from text files
 */
package lapsolver.util;

import lapsolver.EdgeList;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

public class EdgeListLoader {

    public EdgeList loadFromCSV(String fileName) throws IOException {
        return loadFromCSV(fileName, ",");
    }

    public EdgeList loadFromCSV(String fileName, String delimiter) throws IOException {
        ArrayList<Integer> arU = new ArrayList<>();
        ArrayList<Integer> arV = new ArrayList<>();
        ArrayList<Double> arW = new ArrayList<>();

        try (BufferedReader in = new BufferedReader(new FileReader(fileName))) {
            String line = in.readLine();

            while (line != null) {
                String[] values = line.split("\\s*" + Pattern.quote(delimiter) + "\\s*");

                if (values.length != 3)
                    throw new RuntimeException("Error: need exactly three values for edge list entry (u,v,w)");

                arU.add((int) Double.parseDouble(values[0]));
                arV.add((int) Double.parseDouble(values[1]));
                arW.add(Double.parseDouble(values[2]));

                line = in.readLine();
            }
        }

        EdgeList edges = new EdgeList(arU.size());
        for (int i = 0; i < edges.ne; i++) {
            edges.u[i] = arU.get(i);
            edges.v[i] = arV.get(i);
            edges.weight[i] = arW.get(i);
        }
        return edges;
    }
}
