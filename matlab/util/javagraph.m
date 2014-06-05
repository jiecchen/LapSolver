function g = javagraph( a )
%JAVAGRAPH Converts an adjacency matrix to a Java graph.
    import lapsolver.WeightedGraph;
    [ai,aj,av] = find(tril(a));
    g = WeightedGraph();
    g.fromMatlab(ai, aj, av);
end

