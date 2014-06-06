function g = javagraph( a )
%JAVAGRAPH Converts an adjacency matrix to a Java graph.
    import lapsolver.Graph;
    import lapsolver.EdgeList;
    [ai,aj,av] = find(tril(a));
    g = Graph(ai-1, aj-1, av);
end

