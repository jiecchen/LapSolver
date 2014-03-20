function [x] = cgSolver(la,b,opts)
% function [x] = cgSolver(la,b,opts)
% function [f] = cgSolver(la,[],opts)
%
% this calls conjugate gradient with a diagonal preconditioner
%
% opts.tol (default 1e-6)
% opts.maxit (default 1000)
%
% Copyright Daniel Spielman, 2013, Yale University.


default('b',[]);

default('opts','tol',1e-6);
default('opts','maxit',1000);

dd = diag(diag(la));

if isempty(b)
    f = @(b)(internal(la,dd,b,opts));
    x = f;
else
    x = internal(la,dd,b,opts);
end

end % main function


function x = internal(la,dd,b,opts)

  [x,flag] = pcg(la, b, opts.tol, opts.maxit, dd);

end
