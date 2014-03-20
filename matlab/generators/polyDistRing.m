function a = polyDistRing(n,k,e)
% function a = polyDistRing(n,k,e)
%
% on a ring of n vertices, 
% chooses k outgoing edges from every vertex,
% according to distribution 1/k^e
%
% where e is by default 1.
% making e larger concentrates it more
%
% has the ring as a base
%
% example: a = polyDistRing(1000,1,1);
%
% going for the Kleinberg-type ring.
%
% Daniel A. Spielman, Yale University

default('e',1);

mymod = @(x,y)(mod(x,y) + y*(mod(x,y)==0));

ai = [1:n]';
aj = mymod([1:n]'+1,n);

for i = 1:k;
  ai = [ai; [1:n]'];
  add = (1./rand(n,1)).^(1/e);
  add = ceil(add);
  aj = [aj;mymod([1:n]'+add,n)];
end

a = sparse(ai,aj,1,n,n);
a = (a + a' > 0);
a = a - diag(diag(a));

  
