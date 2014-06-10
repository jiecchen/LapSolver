function drawsubgraph(a, s)
%DRAWSUBGRAPH Draws a subgraph.
    xy = specxy(a);
    hold on;
    gplot(a, xy, 'k');
    b = a;
    b(s,s) = 0;
    gplot(a-b, xy, 'r');
    hold off;
end
