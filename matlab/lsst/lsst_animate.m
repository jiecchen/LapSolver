function lsst_animate(a, xy, strat, fps)
%LSST_ANIMATE Shows an LSST building animation.
    import lapsolver.lsst.*;
    import lapsolver.util.*;
    
    if nargin < 2
        xy = specxy(a);
    end
    if nargin < 3
        strat = SimulPathTree;
    end
    if nargin < 4
        fps = 30;
    end
    
    g = a2g(a);
    GraphUtils.reciprocateWeights(g);
    edges = strat.getTreeEdges(g);
    gplot(a, xy, 'c');
    hold on;
    
    for i = 1:edges.ne
        gplot([0 1;1 0], xy([1+edges.u(i) 1+edges.v(i)],:));
        pause(1/fps);
    end
end

