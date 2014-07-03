function [x, err] = kmp1(a, d, b)
%KMP1 Test for KMP1 solver.
    import lapsolver.lsst.*;
    import lapsolver.solvers.*;

    g = a2g(a);
    la = lap(a) + diag(d);
    kmp = KMPSolver(StarDecompositionTree, ConjugateGradientSolver(1000,1e-4));
    kmp.init(g,d);
    x = kmp.solve(b);
    err = norm(la*x - b)
end

