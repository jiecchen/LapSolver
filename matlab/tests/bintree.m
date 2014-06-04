function [ g, t ] = bintree( n )
%BINTREE Binary tree + complete graph test case for LCA, stretch, etc.
    % given n, produces a pair with 2^n - 1 nodes
    % returns: g = java WeightedGraph for complete graph
    %          t = java Tree for binary tree
    nv = 2^n - 1;
    ag = ones(nv,nv) - eye(nv);
    
    % build complete binary tree
    ai = 1:2^(n-1)-1;
    aj = [ (2*ai) (2*ai+1) ];
    ai = [ ai ai ];
    at = sparse(ai,aj,ones(length(ai),1),nv,nv);
    at = at+at';
    
    % convert to java objects
    g = javagraph(ag);
    t = javagraph(at);
    t = t.treeToTree();
end

