function m = e2tril( e, n )
%E2TRIL Converts an edge list to a lower triangular matrix.
    
    m = sparse(double(e.u+1), double(e.v+1), e.weight, n, n);
end
