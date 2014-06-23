function [ x ] = posdef(a)
%POSDEF Tests Kelner solver for a reduced positive definite system.
    n = length(a);
    la = lap(a);
    la = la + diag(rand(1,n))/10;
    b = randn(n,1);
    
    bExt = [b; -sum(b)];
    s = sum(la);
    
    aExt = a;
    aExt(n+1,1:n) = s;
    aExt(1:n,n+1) = s;
    
    laExt = lap(aExt);
    
    % xExt = lapAmdSolver(laExt, bExt);
    xExt = kelner(aExt, bExt, 1000000);
    x = xExt(1:n) - xExt(n+1);
    norm(laExt*xExt - bExt)
    norm(la*x-b)
end

