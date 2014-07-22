function [ rA ] = mtGraph(A, k, tol, iters)
%   This function implements the Miller Tolliver graph reweighting
%   algorithm.
%   edgetol should be ~1e-3 so that the eigenvector linear systems don't
%   crash
%   tol can be played with, 1e-1 works fine
    if nargin < 4
        iters = 1;
    end
    
    rA = A;
    
    for iter = 1:iters
        fprintf('iter %d\n', iter);
        
        [f, e] = compF(rA, k);
        rA = reweightFA(rA, f, e);
        disp(e);
        
        [~,~,w] = find(rA);
        fprintf('reweighted, minw=%.6f, maxw=%.6f\n', min(w), max(w));
        
        rA = mtNormalize(A, rA, f, e, k, tol);
        
        % fprintf('done, norm(rA)=%.6f, norm(A)=%.6f done!\n\n', mtNormFA(rA, A, k), mtNormFA(A, A, k));
    end
end

