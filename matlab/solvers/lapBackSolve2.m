function x = lapBackSolve2(b,l,u)
% function x = lapBackSolve2(b,l,u)
%
% assuming that [l,u] are the lu factorization of a laplacian,
% this does the backsolve to find x st (lu)x = b,
% but first makes b sum to 0.
%
% probably called like
%  finv = @(b)(lapBackSolve(b,l,u));
%


n = length(b);

b = b - mean(b);
y = l \ b;

x = zeros(n,1);
x(1:n-1) = u(1:n-1,1:n-1) \ y(1:n-1);
