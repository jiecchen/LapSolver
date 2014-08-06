function a = knnGraph(x,k)
% function a = knnGraph(x,k)
%
% x is n-by-d, return the directed graph of the k nearest nbrs for each x.
% so, (i,j) is an edge if j is one of the k closest to i
%
% This is a directed graph.
%
% Copyright Daniel Spielman, 2013, Yale University.


[n,d] = size(x);

k = min(k,n-1);

list = zeros(k,n);

for i = 1:n,
  y = x - ones(n,1)*x(i,:);
  s = sum(y.^2,2);  
  [val, ord] = sort(s);
  if (val(2) > eps) % no repeat point
    list(:,i) = ord(2:(k+1));
  else
    [junk,start] = min(val < eps);
    list(:,i) = ord([start:(start+k-1)]);
  end
end

a = sparse(list(:), kron([1:n]', ones(k,1)), 1, n, n);



    
