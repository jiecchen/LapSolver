function a = gaussRing(n,k,sig)
% function a = gaussGraph(n,k,sig)
%
% on a ring of n vertices, 
% chooses k outgoing edges from every vertex,
% according to a gaussian dist with std sigma
%
% has the ring as a base
%
% example: a = gaussRing(1000,3,10);
%



mymod = @(x,y)(mod(x,y) + y*(mod(x,y)==0));

ai = [1:n]';
aj = mymod([1:n]'+1,n);

for i = 1:k;
  ai = [ai; [1:n]'];
  add = sig*randn(n,1);
  add = sign(add).*ceil(abs(add));
  aj = [aj;mymod([1:n]'+add,n)];
end

a = sparse(ai,aj,1,n,n);
a = (a + a' > 0);

  
