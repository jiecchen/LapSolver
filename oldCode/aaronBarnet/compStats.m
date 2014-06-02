function [phi, ne, nv] = compStats(a,comp)
% function [phi, ne, nv] = compStats(a,comp)
%
% for each cluster in comp, 
%   compute ne = number edges attached, with multiplicity
%           nv = number vertices,
%           phi = edges out / ne

if (min(comp) == 0),
  comp = comp + 1;
end


minc = min(comp);
maxc = max(comp);

nc = maxc - minc + 1;

phi = zeros(nc,1);
ne = zeros(nc,1);
nv = zeros(nc,1);

edg = (a > 0);

for c = minc:maxc,
   ind = find(comp == c);
   
   nv(c) = length(ind);
   


   chi = (comp == c);
   
   ne(c) = chi' * a * chi/2;
   
   wtinside = chi' * a * chi;
   bdrywt = (1-chi)' * a * chi;
   
   phi(c) = bdrywt / (wtinside + bdrywt);
end
