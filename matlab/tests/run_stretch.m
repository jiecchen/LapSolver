function st = run_stretch(a)
%RUN_LCA Tests old and new stretch implementations. Returns the stretch.
    import lapsolver.lsst.SimulPathTree;
    import lapsolver.algorithms.Stretch;
    g = javagraph(a);
    spt = SimulPathTree();
    tic; tr = spt.solve(g); toc
    tic; tr.compTotalStretch(g)
    toc
    tic; s = Stretch(g,tr); st = s.totalStretch()
    toc
end