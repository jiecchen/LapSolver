function [ nL ] = normLap( A )
%NORMLAP Computes the normalized laplacian of a graph
    n = length(A);
    la = lap(A);
    d_sqinv = sparse(1:n, 1:n, 1./sqrt(sum(A)), n, n);
    nL = d_sqinv * la * d_sqinv;

    % symmetrize
    nL = tril(nL, -1);
    nL = nL + nL' + speye(n);
end

