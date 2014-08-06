function a = ringGraph(n,gens)
% function a = ringGraph(n,gens)
%
% the cayley graph on n vertices with generators in gens
% 1 is automatically a generator
%
% example: a = ringGraph(n,ceil(sqrt(n)*rand(1,4)));
%
% Daniel A. Spielman, Yale University

mymod = @(x,y)(mod(x,y) + y*(mod(x,y)==0));

ai = [1:n]';
aj = mymod([1:n]'+1,n);

for i = 1:length(gens);
  ai = [ai; [1:n]'];
  aj = [aj;mymod([1:n]'+gens(i),n)];
end

a = sparse(ai,aj,1,n,n);
a = (a + a' > 0);

  
