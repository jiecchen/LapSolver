function stretch = run_lsst(a)
%RUN_LSST Tests the Java -> MATLAB pipeline.
%   Takes a graph, uses Dan's SimulPathTree
%   to generate a low-stretch spanning tree, then computes its mean
%   stretch. Measures time taken by the edgeGrow step.
%   Sample: time_lsst( del3Graph(10000) )
    import lapsolver.lsst.SimulPathLSST;
    g = javagraph(a);
    spt = SimulPathLSST(g);
    tic; tr = spt.solve(g); toc
    stretch = tr.compTotalStretch(g)/g.ne;
end

