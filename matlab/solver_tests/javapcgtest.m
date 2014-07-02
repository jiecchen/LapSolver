function x = javapcgtest(a, b, d)
%JAVAPCGTEST Test for Java CGSolver.
    import lapsolver.solvers.*;
    
    solver = ConjugateGradientSolver(1000, 1e-2);
    solver.init(a2g(a), d);
    x = solver.solve(b);
end

