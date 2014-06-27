function g = a2g(a, datastruct)
%A2G Converts an adjacency matrix to a Java graph object.
%datastruct can be 'graph', 'tree', or 'edgelist' (default 'graph')
    import lapsolver.*;
    
    if nargin < 2
        datastruct = 'graph';
    end
    
    [ai,aj,av] = find(tril(a));
    edges = EdgeList(ai-1, aj-1, av);
    
    if strcmp(datastruct, 'graph')
        g = Graph(edges);
    elseif strcmp(datastruct, 'tree')
        g = Tree(edges);
    else
        g = edges;
    end
end

