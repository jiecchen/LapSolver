function [ reweightedA ] = reweightFA( A, k, f, e )
%REWEIGHTFA Reweight the given graph (given by A) using the Rfa formula.
%           The eigenvalues are found in v and e.

    n = size(A,1);

    a = zeros(n, n, k);
    b = zeros(n, n, k);
    
    [u,v,w] = find(A);
    m = length(u);
    
    for eigIndex = 1:k
        currU = u;
        currV = v;
        currWa = zeros(1,m);
        currWb = zeros(1,m);
        
        for i = 1:m 
            currWa(i) = 1 / e(eigIndex) * (f(u(i),eigIndex) - f(v(i),eigIndex)) ^ 2;
            currWb(i) = f(u(i),eigIndex) ^ 2 + f(v(i),eigIndex) ^ 2;
        end
        
        a(:,:,eigIndex) = sparse(currU, currV, currWa, n, n);
        b(:,:,eigIndex) = sparse(currU, currV, currWb, n, n);
    end
 
    finU = u;
    finV = v;
    finW = zeros(1,m);
    
    % Reweight the graph with the Rfa
    for i = 1:m
        suma = 0;
        sumb = 0;
        
        for eigIndex = 1:k
            suma = suma + a(u(i), v(i), eigIndex);
            sumb = sumb + b(u(i), v(i), eigIndex);
        end
        
        finW(i) = w(i) * (sumb^2 / (suma^2 + sumb^2));
    end
    
    reweightedA = sparse(finU, finV, finW, n, n);
end

