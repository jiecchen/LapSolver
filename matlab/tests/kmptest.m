function kmptest( a, xy )
%KMPTEST Tests the KMP solver.
    if nargin < 2
        xy = specxy(a);
    end

    import lapsolver.lsst.*;
    import lapsolver.solvers.*;
    
    figure
    gplot(a,xy);
    
    figure
    sdt = StarDecompositionTree;
    kmp = KMPSolver(sdt);
    kmp.init(a2g(a));
    hold on;
    gplot(g2a(kmp.sparsifier),xy,'r');
    gplot(g2a(kmp.spanningTree),xy,'k');
    hold off;
    
    figure
    hold on;
    ar = g2a(kmp.reducedSparsifier);
    nr = length(ar);
    xyp = xy(kmp.reductionPerm(end-nr+1:end)+1,:);
    gplot(ar, xyp);
    hold off;
    
    fprintf('removed %.2f%% of vertices\n', 100 * (1 - kmp.reducedSparsifier.nv / kmp.sparsifier.nv));
end

