function lca = run_lca(a)
%RUN_LCA Tests Tarjan's LCA algorithm. Returns the LCA object.
    init
    [ai,aj,av] = find(tril(a));
    g = WeightedGraph();
    g.fromMatlab(ai, aj, av);
    spt = SimulPathLSST(g);
    tic; tr = spt.edgeGrow; toc
    trt = tr.treeToTree;
    lca = TarjanLCA(trt);
end
