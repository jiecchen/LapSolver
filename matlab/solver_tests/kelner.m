function x = kelner( a, b, iters )
%KELNER Run the Kelner algorithm.
    import lapsolver.lsst.*;
    import lapsolver.solvers.KelnerSolver;
    import lapsolver.algorithms.Stretch;

    g = a2g(a);
    
    strat = StarDecompositionTree;
    solver = KelnerSolver(strat);
    solver.init(g);
    x = solver.solve(b, iters);
end

