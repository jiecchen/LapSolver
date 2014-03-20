function [wt,ind] = minVecCut(a,v)
% function [wt,ind] = sparsecut(a,v)
% 
% given an adj matrix a,
% and a vector v,
%
% order a by v,
% and find the minimum weight cut,
% returning the set in ind.

[vp,perm] = sort(v);
ap = a(perm,perm);

n = length(ap);

cuts = cumsum(sum(tril(ap)-triu(ap)));
cuts = cuts(1:(n-1));

[wt,where] = min(cuts);

ind = perm(1:where);


