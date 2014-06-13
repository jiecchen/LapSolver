function g = a2g( a )
%A2G Converts an adjacency matrix to a Java graph object.
    import lapsolver.Graph;
    import lapsolver.EdgeList;
    [ai,aj,av] = find(tril(a));
    g = Graph(ai-1, aj-1, 1./av);
end

