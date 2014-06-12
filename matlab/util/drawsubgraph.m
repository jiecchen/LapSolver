function drawsubgraph(a, s, color)
%DRAWSUBGRAPH Draws a subgraph.
    xy = specxy(a);
    hold on;
    gplot(a, xy, 'k');
    b = a;
    b(s,s) = 0;
    
    if nargin < 3
        color = 'r';
    end
    gplot(a-b, xy, color);
    scatter(xy(s,1), xy(s,2), color);
    
    hold off;
end
