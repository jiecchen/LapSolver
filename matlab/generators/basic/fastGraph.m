function a = fastGraph(n,type,param)
% function a = fastGraph(n,type,param)
%
% adjacency matrix of
%
% if type == '2',
%   then is approx a 2d mesh
% if type == '3',
%   then is approx a 3d mesh
%


if (nargin < 2),
  type = '2';
end


switch (type),
 case '2',
  m = ceil(sqrt(n));
  n = m^2;

  a = sparse(n,n);
  
  a = a + spdiags(ones(n,1),1,n,n);
  a = a + spdiags(ones(n,1),-1,n,n);
  a = a + spdiags(ones(n,1),m,n,n);
  a = a + spdiags(ones(n,1),-m,n,n);
  
 case '3',
  m = ceil(n^(1/3));
  n = m^3;

  a = sparse(n,n);
  
  a = a + spdiags(ones(n,1),1,n,n);
  a = a + spdiags(ones(n,1),-1,n,n);
  a = a + spdiags(ones(n,1),m,n,n);
  a = a + spdiags(ones(n,1),-m,n,n);
  a = a + spdiags(ones(n,1),m^2,n,n);
  a = a + spdiags(ones(n,1),-m^2,n,n);
  
end

  

