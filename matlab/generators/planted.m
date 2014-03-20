function a = planted(n, d1, d2, multi)
% function a = planted(n, d1, d2, multi)
% function a = planted(n, D, [], multi)
%
% return the adj matrix of a graph on 2n vertices,
% with random graphs of degree d1 in each of two
% components,
% joined by a graph of degree d2 between.
%
% if D is a k-by-k matrix, then create a k-part graph with
% degree D(i,j) between parts i and j.
% only uses the upper triangular part of D.
%
% keep multi edges if multi == 1 (the default)
%
% should have an obvious partitioning if d2 < d1 by much at all.
%
% for example, try n = 1000, d1 = 10, d2 = 5;
%
% Copyright Daniel Spielman, 2013, Yale University.

default('multi',1);

if (prod(size(d1)) == 1)

    a11 = randRegular(n,d1);
    a22 = randRegular(n,d1);

    a12 = sparse(n,n);

    for i = 1:d2,
        p = randperm(n);
        a12 = a12 + sparse(1:n,p,1);
    end

    a = [a11, a12; a12', a22];
else
    D = d1;
    D = triu(D) + triu(D,1)';
    k = length(D);
    a = sparse(k*n,k*n);
    
    % first, off-diag blocks
    for i = 1:k,
        for j = (i+1):k,
            a0 = sparse(n,n);
            for ii = 1:D(i,j)
                p = randperm(n);
                a0 = a0 + sparse(1:n,p,1);
            end
            si = (1:n)+(i-1)*n;
            sj = (1:n)+(j-1)*n;
            a(si,sj) = a0;
        end
    end
    
    a = a + a';
    
    % now, diag blocks
    for i = 1:k,
        a0 = randRegular(n,D(i,i));
        si = (1:n)+(i-1)*n;
        a(si,si) = a0;
    end
    
end

if ~multi
  a = double(a > 0);
end

a = a - diag(diag(a));
