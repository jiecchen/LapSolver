function [x, errs] = kelnerwatch( a, b, iters )
%KELNER Watch iterations of the Kelner algorithm.
%Sample:
%errs = kelnerwatch(grid2(100,100),randb(10000),10000);
%semilogy(errs)
    import lapsolver.lsst.SimulPathTree;
    import lapsolver.solvers.KelnerSolver;
    import lapsolver.algorithms.Stretch;

    g = a2g(a);
    
    solver = KelnerSolver(SimulPathTree);
    solver.init(g);
    solver.solve_init(b);

    la = lap(a);
    errs = zeros(1,iters);
    for i = 1:iters
        solver.solve_iter;
        x = solver.solve_return;
        errs(i) = norm(la*x - b) / norm(b);
    end
end

