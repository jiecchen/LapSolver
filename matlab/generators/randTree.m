function a = randTree( n, excess )
%RANDTREE Generates a random tree on n vertices.
%   excess: the number of off-tree edges to add as well.
    if nargin < 2
        excess = 0;
    end
    
    ai = 2:n;
    aj = zeros(1,n-1);
    ak = ones(1,n-1);

    % construct tree edges
    for i = 1:n-1
        aj(i) = randi(i);
    end
    
    a = sparse(ai,aj,ak,n,n);
    a = a + a';
    
    % construct off-tree edges
    for i = 1:excess
        x = -1; y = -1;
        while x == y || a(x,y) ~= 0
            x = randi(n); y = randi(n);
        end
        a(x,y) = 1;
        a(y,x) = 1;
    end
    
end

