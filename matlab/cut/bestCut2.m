function [ind,val] = bestCut2(a,v,kmin,kmax)
% function [ind,val] = bestCut2(a,v,kmin,kmax)
% 
% given an adj matrix a and a vector v,
% order the vertices by v and return
% the prefix in this order that yields the sparsest cut.
%
% only look at sets of size between kmin and kmax

n = length(a);
    
if (nargin < 3),
    kmin = 1;
end

if (nargin < 4),
    kmax = n;
end


% permute a and v according to this order
[vp,perm] = sort(v,'Descend');
ap = a(perm,perm);

% degrees in this permuted order
pdegs = sum(ap);
totdeg = sum(pdegs);

sumwt = zeros(1,n);
edgesin = zeros(1,n);

apu = triu(ap); % upper triangular part
pdegsin = sum(apu); % in-degrees in prefixes

sumdegs = cumsum(pdegs);
edgesin = cumsum(pdegsin(1:n));
bdry = sumdegs - 2*edgesin;

conducs = bdry ./ min(sumdegs, totdeg-sumdegs);

[val,ind] = min(conducs(kmin:kmax));

if (vp(ind+kmin-1)==0)
    S = [1]; % to prevent returning everything
    val = conducs(1);
else
    S = find(v >= vp(ind+kmin-1));
end

ind = zeros(n,1);
ind(S) = 1;
