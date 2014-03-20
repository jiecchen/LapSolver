function c = sparsity(a,s)
% function c = sparsity(a,s)
% function c = sparsity(a,idx)
%
% c is the sparsity of set s in graph a
% or, of idx vector idx.

n = length(a);
if (length(s) < n),
    idx = ones(length(a),1);
    idx(s) = 0;
    sc = find(idx);
else
    idx = s;
    s = find(1-idx);
    sc = find(idx);
end

c = nnz(a(s,sc))/min(nnz(s),nnz(sc));

