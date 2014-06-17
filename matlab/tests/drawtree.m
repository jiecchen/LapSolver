function [t,xy] = drawtree(a, strat, xyin)
%DRAWTREE Takes a tree strategy, and visualizes the result.
    import lapsolver.lsst.*;
    import lapsolver.algorithms.Stretch;
    
    if nargin < 3
        xy = specxy(a);
    else
        xy = xyin;
    end
    
    if nargin < 2
        strat = StarDecompositionTree;
    end
    
    g = a2g(a);
    tic; t = strat.getTree(g); toc;
    stres = Stretch.compute(g,t);
    stretch = stres.total / g.ne
    
    gplot(g2a(t),xy);
end

