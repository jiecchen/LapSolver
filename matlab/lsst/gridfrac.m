function st = gridfrac( n )
%GRIDFRAC Finds the stretch of the X-shaped fractal on a grid graph.
%   Makes a (grid graph, tree) pair with (2^n)^2 vertices
    import lapsolver.algorithms.Stretch;
    import lapsolver.generators.FractalGridGraph;
    import lapsolver.lsst.SimulPathTree;
    import lapsolver.util.GraphUtils;

    g = a2g( grid2(2^n, 2^n) );
    
    % compute SimulPathTree and stretch
    spt = SimulPathTree;
    t = spt.getTree(g);
    Stretch.compute(g,t).total / g.ne
    
    % compute fractal and stretch
    frac = FractalGridGraph(n);
    ft = GraphUtils.toTree( frac.getFractalTree() );
    st = Stretch.compute(g,ft).total / g.ne;
end

