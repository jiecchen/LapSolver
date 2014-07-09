function [ A ] = sanitize( A, tol )
%SANITIZE Removes matrix entries smaller than tol

    N = size(A, 1);
    for i = 1:N
        for j = 1:N
            if A(i, j) < tol
                A(i, j) = 0;
            end
        end
    end
    
end

