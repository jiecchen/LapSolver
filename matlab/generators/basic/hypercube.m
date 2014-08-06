function a = hypercube(d)
% function a = hypercube(d)
%
% adj matrix of a d dim hypercube
%

a = sparse([0 1; 1 0]);

for i = 1:(d-1),
  k = 2^i;
  D = diag(sparse(ones(1,k)));
  a = [a, D; D, a];
end

