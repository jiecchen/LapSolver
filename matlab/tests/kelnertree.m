function t = kelnertree( n )
%KELNERTREE Test the Kelner tree structure.
    import lapsolver.util.GraphUtils;
    import lapsolver.solvers.kelner.KelnerFlowTree;
    import lapsolver.solvers.kelner.DirectFlowTree;
    import lapsolver.EdgeList;
    
    t = GraphUtils.toTree( a2g(cbt(n)) );
    edges = EdgeList([1,3,1,3],[2,6,5,4],[1,1,1,1]);
    kt = KelnerFlowTree(t, edges);
    dt = DirectFlowTree(t, edges);
    
    kt.setTreeFlows([0 1 2 3 4 5 6 7 8 9 10]);
    dt.setTreeFlows([0 1 2 3 4 5 6 7 8 9 10]);
    
    kt.update(1,1);
    kt.update(0,1);
    kt.update(2,-3);
    kt.update(3,-5);
    [kt.query(0) kt.query(1) kt.query(2) kt.query(3)]
    kf = kt.getTreeFlows;
    
    dt.update(1,1);
    dt.update(0,1);
    dt.update(2,-3);
    dt.update(3,-5);
    [dt.query(0) dt.query(1) dt.query(2) dt.query(3)]
    df = dt.getTreeFlows;
    
    kt.dump
    
    [kf df]
end

