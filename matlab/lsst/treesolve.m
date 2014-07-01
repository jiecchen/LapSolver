function err = treesolve( n )
%TREESOLVE Test the tree solver on the complete binary tree.
%   Input: number of nodes in tree. Returns the norm of the error.
    import lapsolver.solvers.TreeSolver;

    a = cbt(n);
    a(1,2) = 2;
    a(2,1) = 2;
    l = lap(a);
    full(l)
    
    solver = TreeSolver;
    solver.init(a2g(a));
    
    b = randb(n);
    x = solver.solve(b);
    y = pinv(full(l)) * b;
    
    err = norm( l*x - b );
end

