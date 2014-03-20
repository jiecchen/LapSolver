function [a,xy,bdry] = imgGraph2(imgname, n, sig);
% function [a,xy,bdry] = imgGraph2(imgname, n, sig);
%
% create a mesh with approx n points from an image
% put more points near dark areas


if (nargin < 3),
  sig = 1;
end

img = double(imread(imgname));

size(img)

dark = 255-img(:,:,1);
dark = dark / max(dark(:));

[i,j,v] = find(dark);

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

%% taken from delTriangle

name = tempname;
nodename = [name , '.node'];

fp = fopen(nodename,'w');

fprintf(fp,'%d 2 0 0\n',size(pts,1));

for i = 1:n
  fprintf(fp,'%d %f %f \n',i,x(i),y(i));
end

fclose(fp);

command = ['triangle/triangle -c ' , name];

outname = [name , '.1.ele'];

unix(command);

%% now, read the output

display('reading output');


fp = fopen(outname,'r');
[fs] = fscanf(fp,'%d %d %d',3);
nt = fs(1);
tri = zeros(nt,3);
for i = 1:nt,
  [fs] = fscanf(fp,'%d %d %d %d',4);
  tri(i,:) = fs(2:4)';
end

fclose(fp);

polyname = [name , '.1.poly'];
fp = fopen(polyname,'r');
[jnk] = fscanf(fp,'%d %d %d %d',4);
[fs] = fscanf(fp,'%d %d',2);
nbdry = fs(1);
bdry = zeros(1,nbdry);
for i = 1:nbdry,
  [fs] = fscanf(fp,'%d %d %d %d',4);
  bdry(i) = fs(2);
end

fclose(fp);


tmp = bdry;
bdry = zeros(1,n);
bdry(tmp) = 1;


display('computing triangles');

a = sparse(n,n);

ii = [tri(:,1); tri(:,2); tri(:,1)];
jj = [tri(:,2); tri(:,3); tri(:,3)];
a = sparse(ii,jj,1,n,n);

a = (a + a') > 0;



a = double(a);

xy = [x(:),y(:)];



delete(nodename);
delete(outname);
delete(polyname);
