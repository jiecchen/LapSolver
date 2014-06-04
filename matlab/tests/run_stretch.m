function st = run_stretch(a)
%RUN_LCA Tests old and new stretch implementations. Returns the stretch.
    import lapsolver.lsst.SimulPathLSST;
    import lapsolver.algorithms.Stretch;
    g = javagraph(a);
    spt = SimulPathLSST(g);
    tic; tr = spt.edgeGrow().treeToTree(); toc
    tic; tr.compTotalStretch(g)
    toc
    tic; s = Stretch(g,tr); st=s.totalStretch()
    toc
end
