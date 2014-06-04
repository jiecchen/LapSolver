function lca = run_lca(a)
%RUN_LCA Tests Tarjan's LCA algorithm. Returns the LCA object.
    import lapsolver.lsst.SimulPathLSST;
    g = javagraph(a);
    spt = SimulPathLSST(g);
    tic; tr = spt.edgeGrow; toc
    trt = tr.treeToTree;
    lca = TarjanLCA(trt);
end
