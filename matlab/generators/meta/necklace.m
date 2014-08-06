function aneck = necklace(a, k, j)
% function aneck = necklace(a, k, j)
%
% take graph a and make a necklace with k copies of a,
% each connected to the next by j random edges.
%
% default for j is 1 and for k is 3
%
% Copyright Daniel Spielman, 2013, Yale University.

default('j',1);
default('k',3);


n = length(a);

aneck = sparse(n*k,n*k);

for i = 1:k,
    p = randperm(n);
    s = (1:n)+(n*(i-1));
    aneck(s,s) = a(p,p);
end

for i = 1:(k-1),
    p1 = randperm(n);
    p2 = randperm(n);
    s1 = p1(1:j)+(n*(i-1));
    s2 = p2(1:j)+(n*(i));
    aneck(s1,s2) = eye(j);
    aneck(s2,s1) = eye(j);
end


p1 = randperm(n);
p2 = randperm(n);
s1 = p1(1:j);
s2 = p2(1:j)+(n*(k-1));
aneck(s1,s2) = eye(j);
aneck(s2,s1) = eye(j);
