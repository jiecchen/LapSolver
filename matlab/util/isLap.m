function v = isLap(m)
% function v = isLap(m)
%
% return true if is symmetric and laplacian
% tests for
%  symmetric
%  non-pos off diags
%  row-sums near zero
%



n = length(m);

v = (max(max(abs(m-m'))) < 2*eps);

if (max(triu(m,1)) > 0), v = false; end

if (max(abs(sum(m))) > n*eps), v = false; end



