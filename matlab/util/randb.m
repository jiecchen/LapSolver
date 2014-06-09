function b = randb( n )
%RANDB Returns a random n-vector orthogonal to the all-ones vector.
    b = rand(n,1);
    b = b - mean(b);
end

