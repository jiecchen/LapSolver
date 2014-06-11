function xy = specxy( a )
%SPECXY Gets a spectral embedding of a graph, given its adjacency matrix.
%   Useful for visualizations.
    la = lap(a);
    [v,~] = eigs(la,3,'sa');
    xy = v(:,2:3);
end

