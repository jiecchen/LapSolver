function g2csv(a, filename, delim)
%G2CSV Exports a graph or adjacency matrix to CSV.
%Accepts a Graph, Tree, EdgeList, or adjacency matrix.
    if isjava(a)
        a = g2a(a);
    end
    if nargin < 3
        delim = ',';
    end
    [ai,aj,av] = find(tril(a,-1));
    dlmwrite(filename,[ai-1 aj-1 1./av], 'delimiter', delim, 'precision', 7);
end
