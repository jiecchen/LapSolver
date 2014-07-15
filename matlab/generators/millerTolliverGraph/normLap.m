function [ nL ] = normLap( A )
%NORMLAP Computes the normalized laplacian of a graph

    n = size(A,1);
        
    % Compute the vertex degrees
    deg = zeros(1, n);
    
    [u, v, w] = find(A);
    m = length(u);
    
    for i = 1:m
        deg(u(i)) = deg(u(i)) + w(i);
        
    finU = u;
    finV = v;
    finW = zeros(1,m);
        
    for i = 1:m
        if u(i) == v(i)
            finW(i) = 1;
        else
            finW(i) = -w(i) / sqrt(deg(u(i)) * deg(v(i)));
        end
    end
    
    nL = sparse(finU, finV, finW, n, n);
end

