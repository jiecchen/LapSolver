function cm = countMatrix(n, h)
% function cm = countMatrix(n, h)
% create a sparse n by h count matrix

i = 1:n;
j = randi(h, 1, n);

v = randi(2, 1, n) * 2 - 3;
cm = sparse(i, j, v, n, h);
end
