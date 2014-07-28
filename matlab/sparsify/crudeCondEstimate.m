function kappa = condEstimate(a,b,solvea,solveb)
% function kappa = condEstimate(a,b,solvea,solveb)
%
% estimate the relative condition number between a and b,
% using solvers for a and b.
%
% a and b should be positive definite
% if no solvers are given, default to incomplete cholesky preconditioners
%


if (isLap(a)),
  defaultStr('solvea','lapIccSolver(a)');
else
  defaultStr('solvea','iccSolver(a)');
end

if (isLap(b)),
  defaultStr('solveb','lapIccSolver(b)');
else
  defaultStr('solveb','iccSolver(b)');
end

n = length(a);

g = @(x)(real(a*solveb(x)));

x = randn(n,1);
for i = 1:100,
    x = x / norm(x);
    x = g(x);
end
e1 = norm(x)




g = @(x)(real(b*solvea(x)));

x = randn(n,1);
for i = 1:100,
    x = x / norm(x);
    x = g(x);
end
e2 = norm(x)


kappa = e1*e2;

