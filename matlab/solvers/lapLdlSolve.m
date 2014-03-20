function x = lapLdlSolve(b,l,d)
% function x = lapLdlSolve(b,l,d)
% 
% given ldl factorization, applies it to b
% modeled on lup solve, but it is missing a perm.
%
% for the Laplacian case

n = length(d);

b = b - mean(b);

y = l \ b;

% if last block is 2x2, use stabler routine
if (nnz(d) == length(d)+2)
  d1 = d;
  d1([n-1,n],[n-1,n]) = eye(2);
  y2 = d1 \ y;
  y2([n-1,n],:) = pinv(full(d([n-1,n],[n-1,n])))*y2([n-1,n],:);
else
  y2 = d \ y;
end

y3 = l' \ y2;

x = y3;

x = x - mean(x);
