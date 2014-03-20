function a = cayleyTorus(n,gens)
% function a = cayleyTorus(n,gens)
%
% the cayley graph on the n-by-n torus,
% with generators in a k-by-2 matrix gens
%
% 
%
% Daniel A. Spielman, Yale University

pair = @(x,n)([floor(x/n), mod(x,n)]);
unpair =  @(xy,n)(n*xy(:,1) + xy(:,2));


ai = [];
aj = [];

k = size(gens,1);

verts = [0:(n^2-1)]';

for i = 1:k;
    ai = [ai; verts];
    thisaj = pair(verts,n) + ones(n^2,1)*gens(i,:);
    thisaj = unpair(mod(thisaj,n),n);
    aj = [aj; thisaj];
end

a = sparse(ai+1,aj+1,1,n^2,n^2);
a = (a + a' > 0);

  
