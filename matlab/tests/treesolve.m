function err = treesolve( n )
%TREESOLVE Test the tree solver on the complete binary tree.
%   Input: number of nodes in tree. Returns the norm of the error.
    import lapsolver.solvers.TreeSolver;

    a = cbt(n);
    l = lap(a);
    
    solver = TreeSolver;
    solver.init(a2g(a));
    
    b = rand(n,1);
    b = b - mean(b)
    x = solver.solve(b);
    
    x = x - x(1)
    y = pinv(full(l)) * b;
    y = y - y(1)
    
    err = norm( l*x - b );
end

