function [a,xy] = delGraph(n, type, param)
% function [a,xy] = delGraph(n, type, param)
%
% type is 'plain', 'gauss', 'skew', 'aniso' 
%
% example a = delGraph(10000,'plain');
%
% or, to see what you get
%
% [a,xy] = delGraph(100,'plain');
% gplot(a,xy)
%


if (nargin < 2)
  type = 'plain';
end

switch (type),
 case 'plain',
  x = rand(1,n);
  y = rand(1,n);
 case 'gauss',
  x = randn(1,n);
  y = randn(1,n);
 case 'skew',
  x = rand(1,n);
  y = rand(1,n)*param;
 case 'aniso',
  x = rand(1,n);
  y = rand(1,n);
end

tri = delaunay(x,y);
a = sparse(n,n);

nt = size(tri,1);

if (max(type == 'aniso') == 1)
  for i = 1:nt,
    a(tri(i,1),tri(i,2)) = cos(atan((x(tri(i,1))-x(tri(i,2))) / ...
				    (y(tri(i,1))-y(tri(i,2)))));
    a(tri(i,1),tri(i,3)) = cos(atan((x(tri(i,1))-x(tri(i,3))) / ...
				    (y(tri(i,1))-y(tri(i,3)))));
    a(tri(i,3),tri(i,2)) = cos(atan((x(tri(i,3))-x(tri(i,2))) / ...
				    (y(tri(i,3))-y(tri(i,2)))));
  end
  a = a + a';
else
  a = sparse(tri(:,1),tri(:,2),1,n,n) + ...
      sparse(tri(:,1),tri(:,3),1,n,n) + ...
      sparse(tri(:,2),tri(:,3),1,n,n);
  a = (a + a') > 0;
end


a = double(a);

xy = [x;y]';
