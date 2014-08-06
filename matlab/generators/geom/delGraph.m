function [a, xy] = delGraph(n, xy)

if nargin < 2
  xy = randxy(n);
end

tri = delaunay(xy(:,1), xy(:,2));

a12 = sparse(tri(:,1),tri(:,2),1,n,n);
a13 = sparse(tri(:,1),tri(:,3),1,n,n);
a23 = sparse(tri(:,2),tri(:,3),1,n,n);

a = double(a12 | a13 | a23);
a = a + a';