function [ mtnorm ] = mtNormFA( mat, initMat, k )
%MTNORM Computes the norm of a matrix for the FA scheme

    % Compute the f vectors for new matrix
    [f, e] = compF(mat, k);
    
    % Get the initial matrix edges
    [u, v, w] = find(initMat);
    m = size(u);
    
    topsum = 0;
    botsum = 0;
    
    for eigIndex = 1:k
        for i = 1:m
            topsum = topsum + w(i) * (f(u(i), eigIndex) - f(v(i), eigIndex)) ^ 2;
            botsum = botsum + w(i) * (f(u(i), eigIndex) ^ 2 + f(v(i), eigIndex) ^ 2);
        end
    end
    
    mtnorm = topsum / botsum;
end

