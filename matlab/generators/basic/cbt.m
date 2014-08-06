function a = cbt(n)
% function a = cbt(n)
%
% make a complete binary tree with n vertices
%

k = floor((n-1)/2);

a = sparse(n,n);

for i = 1:k,
    a(i,2*i) = 1;
    a(i,2*i+1) = 1;
end

if ((2*k+1) < n)
    a(n-1,n) = 1;
end

a(n,n) = 0;
a = a + a';
