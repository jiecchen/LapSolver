function Z = approxReff(a,solver,c)
% function Z = approxReff(a,solver,c)
%
% return a set of vectors Z whose difference norms
% approximate Reffs.
%
% solver should be a solver for the laplacian of a.
%
% reffs should be norm(Z(i,:)-Z(j,:))^2
%



default('c',10);

defaultStr('solver','lapCmgSolver(lap(a))');


[ai,aj,av] = find(triu(a));
Wh = diag(sparse(av.^(1/2)));
U = e2m([ai,aj],-1);
m = length(ai);

nvecs = ceil(c*log(m));

R = randn(m,nvecs);
B = U*Wh*R;

Z = zeros(size(B));

for i = 1:nvecs,
    Z(:,i) = solver(B(:,i));
end

Z = Z / sqrt(nvecs);
