function [f,e] = compF( A, k )
%COMPG Computes the generalized eigenvectors for the pair (L,D) of matrix A
    n = size(A, 1);
    
    % Compute the normalized Laplacian
    nL = normLap(A);
    
    % Compute k eigenvectors for the normalized Laplacian
    
    [allV, allE] = eigs(nL, k + 1, 'sa');
    allE = diag(allE);
   
    e = zeros(k, 1);
    g = zeros(size(nL, 1), k);
    for i = 1:k
        e(i) = allE(i + 1);
        g(:,i) = allV(:,i);
    end
    
    % For each i from 1 to k, compute the f(i) based on g(i)
    d = diag(lap(A));
    for i = 1:n
        d(i) = 1 / d(i) ^ 0.5;
    sqrtInvD = diag(d);
    
    % The f values will solve the system D^(1/2) * f = g.
    % f = D^(-1/2) * g
    
    f = g;
    for i = 1:k
       f(:,i) = sqrtInvD * g(:,i);
    end
end

