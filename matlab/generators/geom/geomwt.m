function [ b ] = geomwt( a, xy, p )
%GEOMWGT re-weight a by the p-norm of its vertices in the xy embedding
%   p=2 by default
if nargin < 3
    p = 2;
end

[i, j] = find(a);
norms = sum(abs(xy(i,:) - xy(j,:)) .^ p, 2) .^ (-1/p);
n = length(a);
b = sparse(i, j, norms, n, n);

end
