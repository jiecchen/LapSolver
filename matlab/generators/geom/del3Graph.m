function [a,xyz] = del3Graph(n, xyz)
% function [a,xyz] = del3Graph(n, type, param)
%
% type is 'plain', 'gauss', 'skew', 'aniso' 
%
% Compute the Delaunay Graph of random points in 3 dimensions, and
% return those points
%

if (nargin < 2)
  xyz = randxy(n, 'uniform', 3);
end

tri = delaunay(xyz(:,1), xyz(:,2), xyz(:,3));

a12 = sparse(tri(:,1),tri(:,2),1,n,n);
a13 = sparse(tri(:,1),tri(:,3),1,n,n);
a14 = sparse(tri(:,1),tri(:,4),1,n,n);
a23 = sparse(tri(:,2),tri(:,3),1,n,n);
a24 = sparse(tri(:,2),tri(:,4),1,n,n);
a34 = sparse(tri(:,3),tri(:,4),1,n,n);

a = double(a12 | a13 | a14 | a23 | a24 | a34);
a = a + a';
