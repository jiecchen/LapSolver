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
e1 = eigs(g,n);
e1max = max(e1);

g = @(x)(real(b*solvea(x)));
e2 = eigs(g,n);
e2max = max(e2);

kappa = e1max*e2max;

