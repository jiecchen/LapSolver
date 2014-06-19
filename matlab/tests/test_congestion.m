function tc = test_congestion(tree)
%   Test if the congestion implementaion is correct

    
    g = javagraph(tree);
    
    
    spt = SimulPathLSST(g);
    tic; tr = spt.edgeGrow().treeToTree(); toc
    tic; tr.compTotalStretch(g)
    toc
    tic; s = Stretch(g,tr); st=s.totalStretch()
    toc
    
end
