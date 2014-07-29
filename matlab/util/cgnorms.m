function norms = cgnorms( solver )
%WATCHNORMS Grab the array of norms from a ConjugateGradientSolver.
    norms = solver.watchNorms(1:solver.watchIters);
end
