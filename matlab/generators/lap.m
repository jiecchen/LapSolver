function l = lap(a)
% function l = lap(a)
%
% given an adjacency matrix a, make the corresponding laplacian.

l = diag(sum(a)) - a;
