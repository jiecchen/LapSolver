function times = kelnertime(a, b, iters, trials)
%KELNERTIME Measure time performance of the Kelner solver.
    import lapsolver.solvers.*;
    import lapsolver.lsst.*;
    
    if nargin < 4
        trials = 10;
    end
    
    ks = KelnerSolver(StarDecompositionTree);
    ks.init(a2g(a));
    
    times = zeros(trials,1);
    for i=1:trials
        tic;
        ks.solve(b,iters);
        times(i) = toc;
    end
end

