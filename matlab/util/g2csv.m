function g2csv(a, filename)
%G2CSV Exports a graph or adjacency matrix to CSV.
%Accepts a Graph, Tree, EdgeList, or adjacency matrix.
    if isjava(a)
        a = g2a(a);
    end
    [ai,aj,av] = find(tril(a,-1));
    csvwrite(filename,[ai-1 aj-1 1./av]);
end
