function stretch = run_lsst(a)
%RUN_LSST Tests the Java -> MATLAB pipeline.
%   Takes a graph, uses Dan's SimulPathTree
%   to generate a low-stretch spanning tree, then computes its mean
%   stretch. Measures time taken by the edgeGrow step.
%   Sample: time_lsst( del3Graph(10000) )
    init
    [ai,aj,av] = find(tril(a));
    g = WeightedGraph();
    g.fromMatlab(ai, aj, av);
    spt = SimulPathLSST(g);
    tic; tr = spt.edgeGrow; toc
    trt = tr.treeToTree;
    stretch = trt.compTotalStretch(g)/length(ai);
end

