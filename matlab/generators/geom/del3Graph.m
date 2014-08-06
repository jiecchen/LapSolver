function [a,xyz] = del3Graph(n, type, param)
% function [a,xyz] = del3Graph(n, type, param)
%
% type is 'plain', 'gauss', 'skew', 'aniso' 
%
% Compute the Delaunay Graph of random points in 3 dimensions, and
% return those points
%


if (nargin < 2)
  type = 'plain';
end

switch (type),
 case 'plain',
  x = rand(1,n);
  y = rand(1,n);
  z = rand(1,n);
 case 'gauss',
  x = randn(1,n);
  y = randn(1,n);
  z = randn(1,n);
 case 'skew',
  x = rand(1,n);
  y = rand(1,n)*param;
  z = rand(1,n);
 case 'aniso',
  x = rand(1,n);
  y = rand(1,n);
  z = rand(1,n);
end

tri = delaunay(x,y,z);
a = sparse(n,n);

nt = size(tri,1);

%if (max(type == 'aniso') == 1)
%  for i = 1:nt,
%    a(tri(i,1),tri(i,2)) = cos(atan((x(tri(i,1))-x(tri(i,2))) / ...
%				    (y(tri(i,1))-y(tri(i,2)))));
%    a(tri(i,1),tri(i,3)) = cos(atan((x(tri(i,1))-x(tri(i,3))) / ...
%				    (y(tri(i,1))-y(tri(i,3)))));
%    a(tri(i,3),tri(i,2)) = cos(atan((x(tri(i,3))-x(tri(i,2))) / ...
%				    (y(tri(i,3))-y(tri(i,2)))));
%  end
%  a = a + a';
%else
  a = sparse(n,n);
  a = a + sparse(tri(:,1),tri(:,2),ones(nt,1),n,n);
  a = a + sparse(tri(:,1),tri(:,3),ones(nt,1),n,n);
  a = a + sparse(tri(:,1),tri(:,4),ones(nt,1),n,n);
  a = a + sparse(tri(:,2),tri(:,3),ones(nt,1),n,n);
  a = a + sparse(tri(:,2),tri(:,4),ones(nt,1),n,n);
  a = a + sparse(tri(:,3),tri(:,4),ones(nt,1),n,n);
  a = (a + a') > 0;
%end


a = double(a);

xyz = [x;y;z]';

ind = find(sum(a));
a = a(ind,ind);
xyz = xyz(ind,:);
