function [v, e] = lapEigs(a, k)
%LAPEIGS Gets the first k eigenvectors
% of the Laplacian matrix, given the adjacency matrix.
    import lapsolver.solvers.*;
    import lapsolver.lsst.*;
    import lapsolver.algorithms.*;
    
    lapsolver = SolverWrapper(ConjugateGradientSolver(100,1e-8),true);
    ipe = InversePowerEigensolver(a2g(a),lapsolver,false);
    
    n = length(a);
    v = [ones(n,1)/sqrt(n) ipe.getVectors(k-1,100)'];
    e = diag(sqrt(sum( (lap(a)*v).^2 )));
end

