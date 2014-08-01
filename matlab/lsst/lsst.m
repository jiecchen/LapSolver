function [t, avgst, allst] = lsst(g, strat)
%LSST Gets a low-stretch spanning tree and its average stretch.
% t = the Tree object
% avgst = average stretch
% g = graph or adjacency matrix
% strat = spanning tree strategy, default StarDecompositionTree
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
        allst = st.allStretches;
    end
end

