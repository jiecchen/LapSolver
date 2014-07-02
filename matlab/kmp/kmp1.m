function [x, err] = kmp1(a, b)
%KMP1 Test for KMP1 solver.
    import lapsolver.lsst.*;
    import lapsolver.solvers.*;

    g = a2g(a);
    D = ones(1,length(a));
    la = lap(a) + diag(D);
    kmp = KMPSolver(StarDecompositionTree);
    x = kmp.solve(g,b,5,D);
    err = norm(la*x - b)
end

