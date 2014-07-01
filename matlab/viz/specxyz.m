function xyz = specxyz( a )
%SPECXYZ Gets a 3D spectral embedding of a graph, given its adjacency matrix.
%   Slightly useful for visualizations.
    la = lap(a);
    [v,~] = eigs(la,4,'sa');
    xyz = v(:,2:4);
end

