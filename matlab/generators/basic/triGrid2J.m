function [ a, coords ] = grid2J( m, n, vw )
% function [a,corrds] = grid2J(m, n, vw)
%
% an m x n grid graph
% if n is not specified, n = m;
%
% vw is degree of anisotropy
% vw is 1 by default
%
    import lapsolver.generators.TriangleGrid2;
    
    if (nargin < 3)
        vw = 1;
    end

    if (nargin < 2)
        n = m;
    end
    
    gen = TriangleGrid2(n,m,vw);
    a = g2a(gen.generateGraph);
    
    coords = double([mod(0:m*n-1, m); idivide(int32(0:m*n-1),int32(m))]');
end

