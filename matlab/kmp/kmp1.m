function [x, err, kmp] = kmp1(a, d, b)
%KMP1 Test for KMP1 solver.
    import lapsolver.lsst.*;
    import lapsolver.solvers.*;

    g = a2g(a);
    la = lap(a) + sparse(1:g.nv, 1:g.nv, d);
    kmp = KMPSolver(StarDecompositionTree, ConjugateGradientSolver(10000,1e-4));
    kmp.init(g,d);
    x = kmp.solve(b);
    err = norm(la*x - b)
end

