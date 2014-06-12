function [x, errs] = kelnerwatch( a, b, iters )
%KELNER Watch iterations of the Kelner algorithm.
%Sample:
%errs = kelnerwatch(grid2(100,100),randb(10000),10000);
%semilogy(errs)
    import lapsolver.lsst.SimulPathTree;
    import lapsolver.solvers.kelner.KelnerSolver;
    import lapsolver.algorithms.Stretch;

    g = a2g(a);
    
    strat = SimulPathTree;
    solver = KelnerSolver(strat);
    solver.init(g);
    solver.solve_init(b);
    
    stretch = Stretch.compute(g, solver.spanningTree).total
    
    errs = zeros(1,iters);
    for i = 1:iters
        solver.solve_iter;
        x = solver.solve_return;
        errs(i) = norm(lap(a)*x - b);
    end
    
    lap(a)*x-b
end

