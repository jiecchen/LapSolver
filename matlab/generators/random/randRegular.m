function a = randRegular(n, k, multi)
% function a = randRegular(n, k, multi)
%
% generate sparse random graph on n nodes of degree approximately k
% 
% if multi = 0, then delete multiedges (default is 1)
%
% Copyright Daniel Spielman, 2013, Yale University.

default('multi',1);

a = sparse(n,n);

for i = 1:k,
  a = a + randMatching(n);
end

a = a + a';

if (multi == 0),
  a = (a + a') > 0;
  a = double(a);
end

  
