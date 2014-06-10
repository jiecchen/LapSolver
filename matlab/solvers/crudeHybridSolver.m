function [x,iter] = crudeHybridSolver(la,perm,k,b,opts)
% function [x,iter] = crudeHybridSolver(la,perm,k,b,opts)
% function [f] = crudeHybridSolver(la,perm,k,[],opts)
%
% This is a crude, two-level hybrid solver.
% It uses a direct outer method and cmg as the inner method.
% Is re-orders the vertices with the order given by perm,
% and then eliminates the first k by LU.
% cmg is used on the rest
%
% It is intended to be used on an ultra-sparsifier (tree + edges)
%
% opts.tol (default 1e-6)
% opts.maxit (default 1000)
%
    
default('b',[]);

laperm = la(perm,perm);

[L,D] = myLDL2(laperm,k);
nnz(L)
Ltrans = L';

Dsub = D(k+1:end,k+1:end);
DsubSolve = cmgSolver(Dsub);

Dtop = diag(D(1:k,1:k));

n = length(la);


if isempty(b)
    f = @(b)(internal(L,Ltrans,Dtop,DsubSolve,perm,k,b));
    x = f;
else
    x = internal(L,Ltrans,Dtop,DsubSolve,perm,k,b);
end

end % main function

function x = internal(L,Ltrans,Dtop,DsubSolve,perm,k,b)
    
  n = length(perm);
  
  bp = b(perm);
  z = L \ bp;
  y = [z(1:k) ./ Dtop; DsubSolve(z(k+1:n))];
  xp = Ltrans \ y;
  x(perm) = xp;
  x = x(:);

end
