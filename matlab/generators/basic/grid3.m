function [a,bdry] = grid3(k1,k2,k3,ani)
% function [a,bdry] = grid3(k1,k2,k3,ani)
%
% a k1 x k2 grid graph
% if k2 is not specified, k2 = k1;
% if k3 is not specified, k3 = k2;
%

if (nargin < 2)
  k2 = k1;
end

if (nargin < 3)
  k3 = k2;
end

if (nargin < 4)
  ani = 1;
end


n = k1*k2*k3;

a = tril(grid2(k1,k2));

a = kron(speye(k3),a);

a = a + ani*spdiags(ones(1,n)', -k1*k2, n, n);

a = a + a';

tm = zeros(k1,k2,k3);
tm(1,:,:) = 1;
tm(end,:,:) = 1;
tm(:,1,:) = 1;
tm(:,end,:) = 1;
tm(:,:,1) = 1;
tm(:,:,end) = 1;
bdry = (tm(:) ==1);
