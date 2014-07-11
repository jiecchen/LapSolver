function [ reweightedA ] = reweightFA( A, k, f, e )
%REWEIGHTFA Reweight the given graph (given by A) using the Rfa formula.
%           The eigenvalues are found in v and e.

    N = size(A,1);

    a = zeros(k, N, N);
    b = zeros(k, N, N);
    
    for eigIndex = 1:k,
        for i = 1:N,
            for j = 1:N,
                if A(i,j) > 0
                    a(eigIndex, i, j) = 1 / e(eigIndex) * (f(i,eigIndex) - f(j,eigIndex)) ^ 2;
                    b(eigIndex, i, j) = f(i,eigIndex) ^ 2 + f(j,eigIndex) ^ 2;
                end
            end
        end
    end
    
    reweightedA = zeros(N);
    
    % Reweight the graph with the Rfa
    for i = 1:N,
        for j = 1:N,
            if A(i,j) > 0
                suma = 0;
                sumb = 0;
                for eigIndex = 1:k,
                    suma = suma + a(eigIndex,i,j);
                    sumb = sumb + b(eigIndex,i,j);
                end
              
                reweightedA(i, j) = A(i,j) * (sumb^2 / (suma^2 + sumb^2));
            end
        end
    end
end

