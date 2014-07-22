function v = normLapEigs(a, k)
%NORMLAPEIGS Gets the first k non-trivial eigenvectors
% of the normalized Laplacian matrix, given the adjacency matrix.
    import lapsolver.solvers.*;
    import lapsolver.lsst.*;
    import lapsolver.algorithms.*;
    
    lapsolver = SolverWrapper(ConjugateGradientSolver(100,1e-8),true);
    nsolver = NormalizedSolver(lapsolver);
    ipe = InversePowerEigensolver(a2g(a),nsolver,true);
    
    v = ipe.getVectors(k,100)';
end

