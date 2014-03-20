function x = ldlSolve(b,l,d)
% function x = ldlSolve(b,l,d)
% 
% given ldl factorization, applies it to b
% modeled on lup solve, but it is missing a perm.

y = l \ b;
y2 = d \ y;
y3 = l' \ y2;

x = y3;
