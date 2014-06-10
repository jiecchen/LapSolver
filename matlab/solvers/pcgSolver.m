function [x,iters] = pcgSolver(la,fpre,b,opts)
% function [x,iters] = pcgSolver(la,fpre,b,opts)
% function f = pcgSolver(la,fpre,[],opts)
%
% This is just a front end for pcg.
% fpre should be a function that approximately solves systems in la
%
% opts.tol (default 1e-6)
% opts.maxit (default 1000)

default('b',[]);

default('opts','tol',1e-6);
default('opts','maxit',1000);

if isempty(b)
    f = @(b)(internal(la,fpre,b,opts));
    x = f;
else
    [x,iters] = internal(la,fpre,b,opts);
end

end % main function

function [x,iters] = internal(la,fpre,b,opts)
    
  [x,flag,relres,iters] = pcg(la, b, opts.tol, opts.maxit,fpre);
    
end
