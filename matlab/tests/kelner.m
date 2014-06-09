function x = kelner( a, b )
%KELNER Summary of this function goes here
%   Detailed explanation goes here
    import lapsolver.lsst.SimulPathTree;
    import lapsolver.solvers.kelner.KelnerSolver;

    strat = SimulPathTree;
    solver = KelnerSolver(strat);
    solver.init(a2g(a));
    
    solver.solve_init(b);
    x = solver.solve_return;
end

