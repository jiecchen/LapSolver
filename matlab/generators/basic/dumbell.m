function a = dumbell(n,p)
% function a = dumbell(n,p)
%
% two n-n-n grids in 3d, joined by a path of length p
%

a = grid3(n,n,n);
[i,j,v] = find(a);

N = length(a);

i2 = i+N;
j2 = j+N;

pi = [1, (2*N+1:2*N+p)]';
pj = [(2*N+1:2*N+p), N+1]';

id = [i;i2;pi];
jd = [j;j2;pj];

a = sparse(id,jd,1,2*N+p,2*N+p);
a = double((a + a') > 0);

