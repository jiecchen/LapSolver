function a = randMatching(n)
% function a = randMatching(n)
%
% return a random matching on n vertices.
% return it as a symmetric matrix
%
% Copyright Daniel Spielman, 2013, Yale University.

p = randperm(n);
n1 = floor(n/2);
n2 = 2*n1;

a = sparse(p(1:n1), p(n1+(1:n1)), 1, n, n);
a = a + a';
