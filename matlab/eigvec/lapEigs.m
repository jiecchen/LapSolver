function v = lapEigs(a, k)
%LAPEIGS Gets the first k non-trivial eigenvectors
% of the Laplacian matrix, given the adjacency matrix.
    import lapsolver.solvers.*;
    import lapsolver.lsst.*;
    import lapsolver.algorithms.*;
    
    lapsolver = SolverWrapper(ConjugateGradientSolver(100,1e-8),true);
    ipe = InversePowerEigensolver(a2g(a),lapsolver,false);
    
    v = ipe.getVectors(k,100)';
end

