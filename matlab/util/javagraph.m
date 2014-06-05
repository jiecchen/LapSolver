function g = javagraph( a )
%JAVAGRAPH Converts an adjacency matrix to a Java graph.
    import lapsolver.Graph;
    [ai,aj,av] = find(tril(a));
    g = Graph();
    g.fromMatlab(ai, aj, av);
end

