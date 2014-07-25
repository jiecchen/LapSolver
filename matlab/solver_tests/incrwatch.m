function [x, errs] = incrwatch(a, b, niters)
%INCRWATCH Watch the incremental meta-solver.
    import lapsolver.solvers.*;
    import lapsolver.lsst.*;
    
    if nargin < 3
        niters = 10; 
    end
    
    % solver = KelnerSolver(StarDecompositionTree);
    solver = KMPSolver(StarDecompositionTree);
    d = ones(length(a),1);
    
    incrsolver = IncrementalSolver(solver, niters);
    incrsolver.init(a2g(a), d);

    incrsolver.solve_init(b);
    errs = zeros(niters,1);
    for i = 1:niters
        res = incrsolver.residue;
        incrsolver.solve_iter;
        x = incrsolver.solve_return;
        errs(i) = norm(incrsolver.residue);
    end

end

