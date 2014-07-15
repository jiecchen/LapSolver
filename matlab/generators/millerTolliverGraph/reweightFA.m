function [ reweightedA ] = reweightFA( A, k, f, e )
%REWEIGHTFA Reweight the given graph (given by A) using the Rfa formula.
%           The eigenvalues are found in v and e.

    n = size(A,1);
    
    [u,v,w] = find(A);
    m = length(u);
    
    finW = zeros(m,1);
    
    for i = 1:m
        suma = 0;
        sumb = 0;
        
        for eigIndex = 1:k
            suma = suma + 1 / e(eigIndex) * (f(u(i),eigIndex) - f(v(i),eigIndex)) ^ 2;
            sumb = sumb + f(u(i),eigIndex) ^ 2 + f(v(i),eigIndex) ^ 2;
        end
        
        finW(i) = w(i) * (sumb^2 / (suma^2 + sumb^2));
    end
    
    reweightedA = sparse(u, v, finW, n, n);
    
    % make reweightedA symmetric
    low = tril(reweightedA, -1);
    reweightedA = low + low';
    
end

