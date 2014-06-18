function m = e2mat( el )
%E2MAT transforms an edge list into a weighted adjacency matrix
%   the size of the matrix (N) is considered as the biggest element in the
%   edgelist

    N = max(max(el.u), max(el.v)) + 1;
    
    A = zeros(N, N);
    
    for i=1:el.ne,
        u = el.u(i:i) + 1;
        v = el.v(i:i) + 1;
        w = el.weight(i:i);
        
        A(u:u,v:v) = A(u:u,v:v) + w;
    end
    
    m = A;
end

