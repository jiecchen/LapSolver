function b = kelner( a, x )
%KELNER Summary of this function goes here
%   Detailed explanation goes here
    import lapsolver.lsst.SimulPathTree;
    import lapsolver.solvers.kelner.KelnerSolver;

    strat = SimulPathTree;
    solver = KelnerSolver(strat);
    solver.init(a2g(a));
    
    solver.solve_init;
    b = solver.solve_return;
end

