function a = erdosRenyi( n, p )
%ERDOSRENYI Generates an Erdos-Renyi random graph G(n,p).
% Hack to ensure connectedness:
    a = tril(rand(n,n) < p,-2);
    a = a + diag(ones(n-1,1),-1);
    a = a + a';
end

