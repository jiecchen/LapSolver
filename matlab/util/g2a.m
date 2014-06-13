function a = g2a( g )
%G2A Converts a Java Graph to an adjacency matrix.
    import lapsolver.EdgeList;
    e = EdgeList(g);
    a = sparse(double(e.u+1), double(e.v+1), 1./e.weight, g.nv, g.nv);
    a = a+a';
end
