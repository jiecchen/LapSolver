function [ reweightedA ] = reweightFA( A, f, e )
%REWEIGHTFA Reweight the given graph (given by A) using the Rfa formula.
%           The eigenvalues are found in v and e.

    n = length(A);
    [u,v,w] = find(A);
    
    sum_a2 = ( (f(u,:) - f(v,:)).^2 * (1./e) ) .^ 2;
    sum_b2 = sum(f(u,:).^2 + f(v,:).^2, 2) .^ 2;
    w_next = w .* (sum_b2 ./ (sum_a2 + sum_b2));
    
    reweightedA = sparse(u, v, w_next, n, n);
    
    % make reweightedA symmetric
    low = tril(reweightedA, -1);
    reweightedA = low + low';
    
end

