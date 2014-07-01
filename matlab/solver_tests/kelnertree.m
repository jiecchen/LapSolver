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
    
    % kt.setTreeFlows(0:(n-1));
    % dt.setTreeFlows(0:(n-1));
    
    %kt.update(1,1);
    %kt.update(0,1);
    kt.update(2,1);
    %kt.update(3,-5);
    kf = kt.getTreeFlows;
    
    %dt.update(1,1);
    %dt.update(0,1);
    dt.update(2,1);
    %dt.update(3,-5);
    df = dt.getTreeFlows;
    
    [kf df]
    nnz(kf-df)
    
    kt.rootStructure.query(2)
end

