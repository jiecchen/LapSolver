function st = run_stretch(a)
%RUN_LCA Tests old and new stretch implementations. Returns the stretch.
    import lapsolver.lsst.SimulPathTree;
    import lapsolver.algorithms.Stretch;
    import lapsolver.algorithms.StretchDan;
    g = a2g(a);
    
    tic; spt = SimulPathTree(g); toc
    tic; st = Stretch.compute(g,spt.getTree).total; toc
end
