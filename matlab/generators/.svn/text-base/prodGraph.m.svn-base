function a = prodGraph(a1,a2)
% function a = prodGraph(a1,a2)
%
% kronecker graph product on adjacency matrices

n1 = length(a1);
n2 = length(a2);
a = kron(speye(n1),a2) + kron(a1,speye(n2));

