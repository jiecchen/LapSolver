function [a,xy,bdry] = imgGraph(imgname, n, sig, exponent);
% function [a,xy,bdry] = imgGraph(imgname, n, sig, exponent);
%
% create a mesh with approx n points from an image
% put more points near high gradient areas
%

default('sig',1);
default('exponent',2);

if (nargin < 3),
  sig = 1;
end

img = double(imread(imgname));

size(img)

dif = img(2:end,2:end,:) - img(1:(end-1),1:(end-1),:);
grad = sum(dif.^exponent,3);
grad = grad / max(grad(:));

[i,j,v] = find(grad);

rat = n/sum(v);

% replace each value in v with a poisson random variable
% to make the expectations correct

numat = random('poiss',v*rat);

% we now make a certain number of points at each position
% we will then perturb them a little bit (by sigma)

pts = [];

ind = find(numat > 0);

while (~isempty(ind))
  
  i = i(ind);
  j = j(ind);
  numat = numat(ind);
  
  pts = [pts;[i,j]];
  
  numat = numat - 1;

  ind = find(numat > 0);
end

pts = pts + randn(size(pts,1),size(pts,2))*sig;

x = pts(:,1);
y = pts(:,2);


n = size(pts,1);

tri = delaunay(x,y);
a = sparse(n,n);

nt = size(tri,1);

  a = sparse(tri(:,1),tri(:,2),1,n,n) + ...
      sparse(tri(:,1),tri(:,3),1,n,n) + ...
      sparse(tri(:,2),tri(:,3),1,n,n);
  a = (a + a') > 0;

  a = double(a);

xy = [x,y];
