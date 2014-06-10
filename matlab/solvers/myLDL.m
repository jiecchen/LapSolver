function [L,D] = myLDL(A,k)
% function [L,D] = myLDL(A,k)
%
% computes L and D so that L*D*L' = A,
% eliminating the first k nodes
%
% for example
%
% >> a = grid2(3);
% >> la = lap(a);
% >> [L,D] = myLDL(la,4);
% >> L = full(L)
%
% it would be interesting to try this with a prefix of the amd or symrcm ordering

n = length(A);

default('k',n);

L = speye(n);
D = sparse(n,n);

for i = 1:k,
  
  m = A((i):end,i) / A(i,i);
  L((i+1):end,i) = m(2:end);
  L(i,i) = 1;
  D(i,i) = A(i,i);

  if (i < n),
    A((i+1):end,i:end) =   A((i+1):end,i:end) - m(2:end)*A(i,i:end);
  end

end

if (k < n);
    D((k+1):n,(k+1):n) = A((k+1):n,(k+1):n);
end


  
