function b = randb( n )
%RANDB Returns a random n-vector orthogonal to the all-ones vector.
    b = rand(1,n);
    b = b - mean(b);
end

