function [ rA ] = mtGraph(A, k, tol, edgetol)
%   This function implements the Miller Tolliver graph reweighting
%   algorithm.
%   edgetol should be ~1e-3 so that the eigenvector linear systems don't
%   crash
%   tol can be played with, 1e-1 works fine
    iter = 0;
    
    while true
        % Compute the f vectors for the A matrix such that 
        % L * f(:,i) = e(i) * D * f(:,i)
        [f, e] = compF(A, k);

        % reweight the matrix
        [u, v, winit] = find(A);
        [u, v, wr] = find(reweightFA(A, f, e));
        wprime = wr;

        n = size(A,1);
        rA = sparse(u, v, wr, n, n);

        alpha = 1;
        while mtNormFA(rA, A, k) > mtNormFA(A, A, k) + tol
            alpha = alpha / 2;
            for i = 1:size(u),
                wprime(i) = (1 - alpha) * winit(i) + alpha * wr(i);
            end

            rA = sparse(u, v, wprime, n, n);
        end

        fprintf('another iteration at level %d\n', iter);
        iter = iter + 1;

        san = sanitize(rA, edgetol);
        if (graphconncomp(san) < k)
            A = rA;
        else
            rA = san;
            break;
        end
    end
end

