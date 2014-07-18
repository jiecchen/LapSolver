function ar = inducedCut(a, colors)
%INDUCEDSUBGRAPHS Separates a graph into its induced subgraphs by color.
    n = length(a);
    [ai,aj,av] = find(a);
    keep = colors(ai) == colors(aj);
    ar = sparse(ai(keep),aj(keep),av(keep),n,n);
end

