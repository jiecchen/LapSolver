function [L,D] = myLDL2(A,k)
% function [L,D] = myLDL2(A,k)
%
% computes L and D so that L*D*L' = A,
% eliminating the first k nodes
%
% This hack leverages lu to work.
% It is absurdly slow, but maybe faster than myLDL.
%
% the result should satsify L*D*L' = A

n = length(A);

default('k',n);

Ahack = A;
Ahack(k+1:n,k+1:n) = 0;
[l,u] = lu(Ahack);
lhack = l;
lhack(1:k,k+1:n) = 0;
lhack(k+1:n,k+1:n) = speye(n-k);
Dhack = lhack \ A / lhack';

L = lhack;
D = Dhack;


  
