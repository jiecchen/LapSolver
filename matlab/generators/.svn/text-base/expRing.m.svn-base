function a = expRing(n,k,sig)
% function a = expGraph(n,k,sig)
%
% on a ring of n vertices, 
% chooses k outgoing edges from every vertex,
% according to an exponential dist with std sigma (approx)
%
% has the ring as a base
%
% example: a = expGraph(1000,3,10);
%


mymod = @(x,y)(mod(x,y) + y*(mod(x,y)==0));

ai = [1:n]';
aj = mymod([1:n]'+1,n);

for i = 1:k;
  ai = [ai; [1:n]'];
  add = -log(1./rand(n,1));
  add = sig*add;
  add = sign(add).*ceil(abs(add));
  aj = [aj;mymod([1:n]'+add,n)];
end

a = sparse(ai,aj,1,n,n);
a = (a + a' > 0);

a = a - diag(diag(a));  
