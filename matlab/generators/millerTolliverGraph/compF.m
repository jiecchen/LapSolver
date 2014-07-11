function [f,e] = compF( A, k )
%COMPG Computes the generalized eigenvectors for the pair (L,D) of matrix A
    N = size(A, 1);
    
    % Compute the normalized Laplacian
    nL = normLap(A);
    
    % Compute k eigenvector for the normalized Laplacian
    % opts.maxit = 1000;
    % [g,e] = eigs(nL, k, 'SM', opts);
    % e = diag(e);
    
%    disp('**********************************************************');
%    graphconncomp(sparse(full(nL)));
%    full(nL);
%    rank(full(nL));
   
    [allV, allE] = eig(full(nL));
    allE = diag(allE);
    
    %disp('*********************** allV and allE');
    %allV
    %allE
   
    e = zeros(k, 1);
    g = zeros(size(nL, 1), k);
    for i = 1:k
        e(i) = allE(i + 1);
        g(:,i) = allV(:,i);
    end
    
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

