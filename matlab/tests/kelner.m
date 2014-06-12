function x = kelner( a, b, iters )
%KELNER Run the Kelner algorithm.
    import lapsolver.lsst.SimulPathTree;
    import lapsolver.solvers.kelner.KelnerSolver;
    import lapsolver.algorithms.Stretch;

    g = a2g(a);
    
    strat = SimulPathTree;
    solver = KelnerSolver(strat);
    solver.init(g);
    x = solver.solve(b, iters);
end

