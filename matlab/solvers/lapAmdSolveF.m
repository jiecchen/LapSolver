function [solvela,l] = lapAmdSolveF(la)
% function [solvela,l] = lapAmdSolveF(la)
%
% la should be a laplacian.
% does a backsolve using lu factorization from amd order,
% and does it by projecting into nullspace
%
% exmple:
% a = grid2(50);
% la = lap(a) + speye(2500)/10^9;
% f = amdSolveF(la);
% b = randn(2500,1);
% norm(la*f(b) - b)
%



p = symamd(la);
plap = la(p,p);
[l,u] = lu(plap);


permp = @(y)(y(p));
invp(p) = [1:length(p)];
ipermp = @(y)(y(invp));

solvela = @(b)(ipermp(lapBackSolve(permp(b),l,u)));
