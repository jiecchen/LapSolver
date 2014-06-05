function [ a, coords ] = grid2J( m, n, vw )
% function [a,corrds] = grid2J(m, n, vw)
%
% an m x n grid graph
% if n is not specified, n = m;
%
% vw is degree of anisotropy
% vw is 1 by default
%
    import lapsolver.generators.Grid2;
    
    if (nargin < 3)
        vw = 1;
    end

    if (nargin < 2)
        n = m;
    end
    
    total = m*n;
    gen = Grid2(m,n,vw);
    ijv = gen.generateGraph.toIJV();
    a = sparse(ijv(1,:)+1, ijv(2,:)+1, ijv(3,:), total, total);
    a = a + a';
    
    coords = double([mod(0:m*n-1, n); idivide(int32(0:m*n-1),int32(n))]');
    coords = coords / max(max(coords));
end

