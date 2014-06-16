function [u, w] = a2u( a )
%A2U Converts an adjacency matrix to a vertex-edge incidence matrix.
    %Each column (representing an edge) has two entries, a 1 and a -1.

    [ai, aj, av] = find(tril(a,-1));
    nv = length(a);
    ne = length(ai);
    
    w = spdiags(av,0,ne,ne);
    u = sparse(ai,1:ne,ones(1,ne),nv,ne) - sparse(aj,1:ne,ones(1,ne),nv,ne);
end

