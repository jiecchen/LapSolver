function [x,iter] = pcgAmdSolver(la,lpre,b,opts)
% function [x] = pcgAmdSolver(la,lpre,b,opts)
% function [f] = pcgAmdSolver(la,lpre,[],opts)
% function [x] = pcgAmdSolver(la,lsst,b,opts)
%
% This calls pcg with preconditioner given by lpre.
% lpre in turn is solved by amdSolver
% For example, lpre could be from a low stretch tree.
% Or, a low stretch tree plus edges
%
% alternatively, it can take a SpanningTreeStrategy, lsst,
% in which case it uses it to construct a tree.
% 
% lpre is required!
%
% opts.tol (default 1e-6)
% opts.maxit (default 1000)
%
    
if isa(lpre,'lapsolver.lsst.SpanningTreeStrategy'),
    lsst = lpre;
    g = a2g(diag(diag(la))-la);
    tr = lsst.getTree(g);
    t = g2a(tr);
    lt = lap(t) + diag(sum(la));
    fpre = amdSolver(lt);
else
    fpre = amdSolver(lpre);
end

    

default('b',[]);

default('opts','tol',1e-6);
default('opts','maxit',1000);


if isempty(b)
    f = @(b)(internal(la,fpre,b,opts));
    [x,iter] = f;
else
    [x,iter] = internal(la,fpre,b,opts);
end

end % main function


function [x,iter] = internal(la,fpre,b,opts)

  [x,flag,relres,iter] = pcg(la, b, opts.tol, opts.maxit, fpre);

end
