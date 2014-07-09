function [f,e] = compF( A, k )
%COMPG Computes the generalized eigenvectors for the pair (L,D) of matrix A

    % Compute the normalized Laplacian
    nL = normLap(A);

    full(nL);
    
    % Compute k eigenvector for the normalized Laplacian
    [g,e] = eigs(nL, k, 'sm');
    e = diag(e);
    
    % For each i from 1 to k, compute the f(i) based on g(i)
    D = zeros(size(A,1));
    for i = 1:size(A,1)
        for j = 1:size(A,1)
            D(i,i) = D(i,i) + A(i,j);
        end
    end
    
    % The f values will solve the system D^(1/2) * f = g.
    sqrtD = D^0.5;
    
    f = g;
    for i = 1:k
        f(:,i) = linsolve(sqrtD, g(:,i));
    end
end

