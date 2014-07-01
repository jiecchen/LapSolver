function drawgraph( a, opt )
%DRAWGRAPH Draws a graph using its spectral embedding.
    if nargin < 2
        gplot(a, specxy(a));
    elseif strcmp(opt,'3d')
        gplot3(a, specxyz(a));
    end
end

