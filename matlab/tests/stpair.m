function [g, t] = stpair( a )
%STPAIR Takes adjacency matrix of a graph, returns graph and spanning tree
%objects
    import lapsolver.lsst.SimulPathTree;
    g = a2g(a);
    spt = SimulPathTree();
    t = spt.solve(g);
end

