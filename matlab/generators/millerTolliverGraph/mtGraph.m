function [ rA ] = mtGraph( A, k, tol, edgetol, level)
%   This function implements the Miller Tolliver graph reweighting
%   algorithm.
%   edgetol should be ~1e-3 so that the eigenvector liniar systems don't
%   crash
%   tol can be played with, 1e-1 works fine

    % Compute the f vectors for the A matrix such that 
    % L * f(:,i) = e(i) * D * f(:,i)
    [f, e] = compF(A, k);

    % reweight the matrix
    [u, v, winit] = find(A);
    [u, v, wr] = find(reweightFA(A, k, f, e));
    wprime = wr;
    
    N = size(A,1);
    rA = sparse(u, v, wr, N, N);
    
    alpha = 1;
    while mtNormFA(rA, A, k) > mtNormFA(A, A, k) + tol
        alpha = alpha / 2;
        for i = 1:size(u),
            wprime(i) = (1 - alpha) * winit(i) + alpha * wr(i);
            disp(wprime(i))
        end
        
        rA = sparse(u, v, wprime, N, N);
    end
    
    fprintf('another iteration at level %d\n', level)
    
    san = sanitize(rA, edgetol);
    if (graphconncomp(san) < k)
        rA = mtGraph(rA, k, tol, edgetol, level + 1);
    else
        rA = san;
    end
end

