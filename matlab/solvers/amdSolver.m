function [x, l, u, p] = amdSolver(la,b,opts)
% function [f, l, u, p] = amdSolver(la)
% function [x, l, u, p] = amdSolver(la,b)
%
% la should be pos def.
% This is a direct solver that uses the amd order.
%
% exmple:
% a = grid2(50);
% la = lap(a) + speye(2500)/10^9;
% f = amdSolver(la);
% b = randn(2500,1);
% norm(la*f(b) - b)
%
% Copyright Daniel Spielman, 2013, Yale University.


default('b',[]);

p = symamd(la);
plap = la(p,p);
[l,u] = lu(plap);

f = lupSolver(l,u,p);

if ~isempty(b)
    x = f(b);
else
    x = f;
end
