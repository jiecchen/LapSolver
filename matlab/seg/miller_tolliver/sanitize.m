function [ A ] = sanitize( A, tol )
%SANITIZE Removes matrix entries smaller than tol

    n = size(A, 1);
    
    [u, v, w] = find(A);
    m = length(u);
    
    for i = 1:m
        if (w(i) < tol)
            w(i) = 0;
        end
    end
    
    A = sparse(u, v, w, n, n);
end

