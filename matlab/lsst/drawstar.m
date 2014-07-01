function drawstar(a, x0, xy)
%DRAWSTAR Draws the star coloring of a graph.
    import lapsolver.lsst.StarDecompositionTree;

    if nargin < 2
        x0 = 0;
    end
    
    if nargin < 3
        xy = specxy(a);
    end
    
    g = a2g(a);
    
    sdt = StarDecompositionTree;
    c = sdt.getStarColoring(g,x0);
    
    hold on;
    gplot(a(c==0,c==0),xy(c==0,:),'k');
    
    colors = ['r', 'y', 'g', 'b', 'm', 'c'];
    
    notcut = nnz(c==-1)
    
    for i = 1:max(c)
        mask = c==i;
        color = colors(1+mod(i,length(colors)));
        gplot(a(mask,mask), xy(mask,:), color);
    end
    hold off;
end

