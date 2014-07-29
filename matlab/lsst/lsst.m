function [t, avgst] = lsst(g, strat)
%LSST Gets a low-stretch spanning tree and its average edge stretch.
    import lapsolver.algorithms.*;
    import lapsolver.lsst.*;
    import lapsolver.util.*;
    
    if nargin < 2
        strat = StarDecompositionTree;
    end
    
    if strcmp(class(g), 'double')
        g = a2g(g);
    end
    
    GraphUtils.reciprocateWeights(g);
    t = strat.getTree(g);
    GraphUtils.reciprocateWeights(g);
    
    if nargout >= 2
        st = Stretch.compute(g,t);
        avgst = st.total / g.ne;
    end
end

