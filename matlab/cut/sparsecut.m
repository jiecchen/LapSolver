function [st,sumwt,edgein,perm] = sparsecut(a,v,k)
% function [st,sumwt,edgein,perm] = sparsecut(a,v,k)
% 
% given an adj matrix a,
% and a vector v,
%
% order a by v,
% and return as a set the sparsest cut found.
%
% do not look at sets of size smaller than k

default('k',1);

[vp,perm] = sort(v);
ap = a(perm,perm);
aps = sum(ap);


n = length(a);

sumwt = zeros(1,n);
edgein = zeros(1,n);

apu = triu(ap);
apus = sum(apu);

%sumwt(1) = sum(ap(1,:));
%for i = 2:n,
%  sumwt(i) = sumwt(i-1) + sum(ap(:,i));
%end
sumwt = cumsum(aps);


%if (norm(sumwt-sumwt2) > eps), error('here'); else display('ok'); end


%edgein = zeros(1,n);
%for i = 2:n,
%  edgein(i) = edgein(i-1) + sum(apu(:,i));
%end
edgein = [0,cumsum(apus(2:n))];

%if (norm(edgein-edgein2) > eps), error('here'); end


nh = n-1-k;

[val,ind] = min((sumwt(k:nh)-2*edgein(k:nh)) ./ (sumwt(k:nh) .* ...
						 (sumwt(end)-sumwt(k:nh))));

st = find(v <= vp(ind+k-1));

