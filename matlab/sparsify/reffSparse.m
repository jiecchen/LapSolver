function lb = reffSparse(la,avedeg,solvela)
% function lb = reffSparse(la,avedeg,solvela)
%
% Given a laplacian matrix la, produce a sparsified version lb,
% by sampling according to effective resistances, following 
% Spielman-Srivastava.
%
% it would like a solver for the laplacian of a.  By default, it uses incomplete
% Cholesky.
%
% avedeg should be at least log n.  It defaults to 8 * log n.
%


n = length(la);

% get the adjacency matrix
at = triu(-la,1);  a = at + at';

default('avedeg',ceil(4*log(n)));

default('solvela',iccSolver(la));

Z = approxReff(a,solvela);

[ai,aj,av] = find(at);

reff = zeros(size(av));
for i = 1:length(reff),
    reff(i) = norm(Z(ai(i),:)-Z(aj(i),:))^2;
end

probs = (reff.*av*avedeg/2);
samp = (rand(size(reff)) < probs);

b = sparse(ai(samp),aj(samp),av(samp)./probs(samp),n,n);
b = b + b';
lb = diag(sum(b)) - b;
