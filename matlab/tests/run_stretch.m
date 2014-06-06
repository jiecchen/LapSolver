function st = run_stretch(a)
%RUN_LCA Tests old and new stretch implementations. Returns the stretch.
    import lapsolver.lsst.SimulPathTree;
    import lapsolver.algorithms.Stretch;
    import lapsolver.algorithms.StretchDan;
    g = a2g(a);
    
    spt = SimulPathTree();
    tic; tree = spt.getTree(g); toc
    tic; st = Stretch.compute(g,tree).total; toc
end
