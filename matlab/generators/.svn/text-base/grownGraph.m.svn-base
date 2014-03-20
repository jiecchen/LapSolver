function a = grownGraph(n,k)
% function a = grownGraph(n,k)
%
% create a graph on n vertices
% for each vertex, give it k edges to randomly chosen prior
% vertices
%

a = sparse(n,n);

for i = 1:k,
  e = [[1:n]', [1;ceil(rand(n-1,1).*([1:(n-1)]'))]];
  a = a + sparse(e(:,1),e(:,2),1,n,n);
end

a = a - diag(diag(a));
a = a + a';



