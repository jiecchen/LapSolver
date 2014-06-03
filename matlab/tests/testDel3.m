function stretch = testDel3()
%TESTDEL3 Tests the Java -> MATLAB pipeline.
%   Generates a 3D Delaunay graph, uses Dan's SimulPathTree
%   to generate a low-stretch spanning tree, then computes its mean
%   stretch.
    a = del3Graph(10000);
    [ai,aj,av] = find(tril(a));
    g = WeightedGraph();
    g.fromMatlab(ai, aj, av);
    spt = SimulPathTree(g);
    tr = spt.edgeGrow;
    trt = tr.treeToTree;
    stretch = trt.compTotalStretch(g)/length(ai);
end

