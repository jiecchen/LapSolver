function [ outG ] = randomizeEdges( graph )
%RANDOMIZEEDGES Given a graph, multiplies the edges by some random quantity

    adjMat = g2a(graph);
    [u,v,w] = find(adjMat);
    for i = 1:length(w)
        w(i) = w(i) * rand(1,1);
    end
    adjMat = sparse(u,v,w);
    
    outG = a2g(adjMat);
end

