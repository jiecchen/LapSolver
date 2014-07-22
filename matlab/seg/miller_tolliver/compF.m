function [f,e] = compF( A, k )
%COMPG Computes the generalized eigenvectors for the pair (L,D) of matrix A
    n = length(A);
    
    % Compute the normalized Laplacian
    nL = normLap(A);
    
    % Compute k eigenvectors for the normalized Laplacian
    
    % opts.tol = 1e-5;
    % [g, e] = eigs(nL, k, 'sa', opts);
    [g, e] = normLapEigs(A, k);
    
    e = diag(e);
    
    % g = g(:,2:(k+1));
    % e = e(2:(k+1));
    
    % The f values will solve the system D^(1/2) * f = g
    % f = D^(-1/2) * g
    f = bsxfun(@times, 1 ./ sqrt(sum(A))', g);
end
