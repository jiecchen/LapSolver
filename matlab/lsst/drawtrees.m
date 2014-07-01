function drawtrees( a, xy )
%DRAWTREES Compares SimulPathTree to StarDecompositionTree
    if nargin < 2
        xy = specxy(a);
    end
    
    import lapsolver.lsst.*;
    
    figure
    disp('StarDecompositionTree')
    drawtree(a, StarDecompositionTree, xy);
    title('StarDecompositionTree')
    
    figure
    disp('SimulPathTree');
    drawtree(a, SimulPathTree, xy);
    title('SimulPathTree')
end

