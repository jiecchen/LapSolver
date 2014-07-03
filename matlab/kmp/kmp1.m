function [x, err] = kmp1(a, b)
%KMP1 Test for KMP1 solver.
    import lapsolver.lsst.*;
    import lapsolver.solvers.*;

    g = a2g(a);
    D = ones(1,length(a));
    la = lap(a) + diag(D);
    kmp = KMPSolver(StarDecompositionTree, ConjugateGradientSolver(1000,1e-4));
    kmp.init(g,D);
    x = kmp.solve(b);
    err = norm(la*x - b)
end

