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

        n = length(A);
        rA = sparse(u, v, wr, n, n);

        fprintf('another iteration at level %d with norms (%.10f, %.10f)\n', iter, mtNormFA(rA, A, k), mtNormFA(A, A, k));
        fprintf('min edge is %.10f and max edge is %.10f\n', min(winit), max(winit));
        
        alpha = 1;
        while mtNormFA(rA, A, k) > mtNormFA(A, A, k) + tol
            alpha = alpha / 2;
            for i = 1:size(u),
                wprime(i) = (1 - alpha) * winit(i) + alpha * wr(i);
            end

            rA = sparse(u, v, wprime, n, n);
        end
        
        fprintf('normalizing done!\n\n');

        iter = iter + 1;
        
        %if (iter == 30)
        %    rA = smartCluster(rA, k);
        %    break;
        %end
        
        if mtNormFA(rA, A, k) < tol
            %sanitize so that only k components remain
            fprintf('**********************************************\n');
            rA = incrementalSanitize(rA, k);
            break;
        end
        
        san = sanitize(rA, edgetol);
        if (graphconncomp(san) < k)
            A = rA;
        else
            rA = san;
            break;
        end
    end
end

