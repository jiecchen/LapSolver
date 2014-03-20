function f = ldlSolveF(l,d)
% function f = ldlSolveF(l,d)
%
% given an ldl factorization, applies it to solve for x
% returns a function f such that f(b) does this

f = @(b)(ldlSolve(b,l,d));
