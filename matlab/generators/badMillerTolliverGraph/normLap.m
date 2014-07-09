function [ nL ] = normLap( A )
%NORMLAP Computes the normalized laplacian of a graph

    N = size(A,1);
        
    % Compute the vertex degrees
    deg = zeros(1, N);

    for i = 1:N,
        for j = 1:N,
            deg(i) = deg(i) + A(i,j);
        end
    end
    
    nL = zeros(N);
    
    for i = 1:N,
        for j = 1:N,
            if i == j
                nL(i,j) = 1;
            elseif A(i,j) > 0
                nL(i,j) = -A(i,j) / sqrt(deg(i) * deg(j));
            end
        end
    end

end

