function [v, e] = normLapEigs(a, k)
%NORMLAPEIGS Gets the first k eigenvectors
% of the normalized Laplacian matrix, given the adjacency matrix.
    import lapsolver.solvers.*;
    import lapsolver.lsst.*;
    import lapsolver.algorithms.*;
    
    lapsolver = SolverWrapper(ConjugateGradientSolver(100,1e-8),true);
    nsolver = NormalizedSolver(lapsolver);
    ipe = InversePowerEigensolver(a2g(a),nsolver,true);
    
    dsqrt = sqrt(full(sum(a))');
    v = [dsqrt/norm(dsqrt) ipe.getVectors(k-1,100)'];
    e = diag(sqrt(sum( (normLap(a)*v).^2 )));
end

