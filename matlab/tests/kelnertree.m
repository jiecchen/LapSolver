function t = kelnertree( n )
%KELNERTREE Test the Kelner tree structure.
    import lapsolver.util.GraphUtils;
    import lapsolver.solvers.kelner.KelnerFlowTree;
    import lapsolver.EdgeList;
    
    t = GraphUtils.toTree( a2g(cbt(n)) );
    edges = EdgeList([1,3],[2,6],[1,1]);
    kt = KelnerFlowTree(t, edges);
    
    kt.update(0,1);
    kt.query(0)
    kt.update(1,-1);
    kt.query(0)

end

