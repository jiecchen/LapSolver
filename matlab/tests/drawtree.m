function [t,xy] = drawtree(a, strat)
%DRAWTREE Takes a tree strategy, and visualizes the result.
    import lapsolver.lsst.*;
    xy = specxy(a);
    if nargin < 2
        strat = StarDecompositionTree
    end
    g = a2g(a);
    t = strat.getTree(g);
    gplot(g2a(t),xy);
end

